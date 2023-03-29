import scipy.io as spio
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
import ionChem
import numpy as np
import matplotlib.pyplot as plt
import loadMSIS
import loadmat
import pickle


def ic(direc, chemistry_config, file, iteration, mixf = 0):
    # load content of last Elspec iteration
    f = direc + file + str(iteration)
    mat = loadmat.loadmat(f)
    con = mat["ElSpecOut"]

    ne = con["ne"]
    n_model = con["iri"]
    [Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
    #normalise charged species to fit electron density
    [nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * ne

    # setting start time [0] to 0:
    ts = con["ts"] - con["ts"][0]
    te = con["te"] - con["ts"][0]
    e_prod = con["q"]
    z_model = con["h"]

    ts_ = np.copy(ts)
    ts[0] = -30*60
    ts_show = np.arange(ts[0], te[-1], 0.01)

    def stepped_prod_t(prod, t):
        """
        returns the production according to the ELSPEC model
        of any species where prod is defined
        at arbitrary times
        using the last ts: ts <= t from ELSPEC
        """
        if t < ts[0]:
            return 0
        if t > te[-1]:
            return 0
        else:
            i_max_ts = len(ts[ts <= t]) - 1
            if i_max_ts < 0:
                print(i_max_ts)
                raise RuntimeError
            prod_t = prod[:, i_max_ts]
            return prod_t

    model = ionChem.ionChem(chemistry_config, z_model)

    for c in model.all_species:
        c.density = nN2 * 0

    # assign densities
    model.e.density = ne
    model.N2.density = nN2
    model.O2.density = nO2
    model.O.density = nO
    model.NOp.density = nNOp
    model.O2p.density = nO2p
    try:
        model.Op.density = nOp
    except AttributeError:
        model.Op_4S.density = nOp

    #custom ion densities:
    #if iteration == 0:
    #    model.NOp.density = model.NOp.density/2
    #    model.O2p.density = model.O2p.density*2
    #
    #    model.O2.density = nO2-nO/2
    #    model.O.density = nO * 2
    #     for n, i in enumerate(model.ions):
    #         i.density[:, 0] = i.density[:, 0] * 0.6
    # model.e.density = np.sum(np.array([i.density for i in model.ions]), axis = 0)


    msis_model = loadMSIS.loadMSIS_new('Data/other/msis.rtf')
    nH = msis_model[4]
    nH_intp = np.exp(PchipInterpolator(msis_model[0][1:-3] / 1e3, np.log(nH[1:-3]))(z_model))
    model.H.density = np.tile(nH_intp, (len(ts), 1)).T

    model.check_chargeNeutrality()

    for c in model.all_species:
        c.prod = e_prod * 0

    Op_prod = e_prod * 0.56 * model.O.density / \
               (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)
    O2p_prod = e_prod * 1.00 * model.O2.density / \
               (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)
    N2p_prod = e_prod * 0.92 * model.N2.density / \
               (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)

    # create smooth function for production
    t = np.arange(0, te[-1], 0.01)

    stepf = np.array([stepped_prod_t(e_prod, i) for i in t])
    e_prod_smooth = PchipInterpolator(t, stepf)

    stepf = np.array([stepped_prod_t(Op_prod, i) for i in t])
    Op_prod_smooth = PchipInterpolator(t, stepf)

    stepf = np.array([stepped_prod_t(O2p_prod, i) for i in t])
    O2p_prod_smooth = PchipInterpolator(t, stepf)

    stepf = np.array([stepped_prod_t(N2p_prod, i) for i in t])
    N2p_prod_smooth = PchipInterpolator(t, stepf)

    # zero function for species that dont feel production
    def zerof(t):
        try:
            some_object_iterator = iter(t)
            return np.zeros(len(t), len(z_model))
        except TypeError:
            return np.zeros(len(z_model))

    # save production function as array, in order of species
    prodMat = np.array([zerof for i in range(len(model.all_species))], dtype='object')

    for i, c in enumerate(model.all_species):
        # print(c.name)
        if c.name == 'e':   prodMat[i] = e_prod_smooth
        if c.name == 'O+':  prodMat[i] = Op_prod_smooth
        if c.name == 'O+(4S)': prodMat[i] = Op_prod_smooth
        if c.name == 'O2+': prodMat[i] = O2p_prod_smooth
        if c.name == 'N2+': prodMat[i] = N2p_prod_smooth

    # new way to do ionospheric chemistry
    # 1. take all info from reaction (reaction rate, constituents)
    ode_mat = np.zeros((len(model.all_reactions), len(model.all_species)), dtype='object')
    rrate = np.array([np.array([r.r_rate_t(*t) for r in model.all_reactions]) for t in zip(Tn.T, Ti.T, Te.T)])

    for r in model.all_reactions:
        ed = r.educts_ID
        pr = r.products_ID
        br = r.branching
        # print(e)
        for ed_i in ed:
            ode_mat[r.r_ID, ed_i] = np.array([r.r_ID, -1, *ed])
        for pr_i, bra in zip(pr, br):
            ode_mat[r.r_ID, pr_i] = np.array([r.r_ID, +1 * bra, *ed])

        for ed_i in ed:
            for pr_i in pr:
                if ed_i == pr_i: ode_mat[r.r_ID, ed_i] = 0

    # 2. produce raw DG from 1., excluding all terms that are 0 anyways.
    ode_raw = np.empty(len(model.all_species), dtype='object')
    for i in model.all_species:
        #    print(np.array([o for o in ode_mat[:, i.c_ID] if type(o)!= int]))
        ode_raw[i.c_ID] = np.array([o for o in ode_mat[:, i.c_ID] if type(o) != int])

    def asint(arr):
        return arr.astype(int)

    # 3. produce DG with only relevant terms
    def fun(t, n, h):
        if ts[0] <= t <= te[-1]:
            k = len(ts[ts <= t]) - 1
        else:
            print(t)
            raise RuntimeError

        dndt = np.array([np.sum(((rrate[k, asint(ode_raw[i.c_ID][:, 0]), h].T * ode_raw[i.c_ID][:, 1]).T \
                                 * n[asint(ode_raw[i.c_ID][:, 2])] * n[asint(ode_raw[i.c_ID][:, 3])] \
                                 ), axis=0) \
                         + prodMat[i.c_ID](t)[h] \
                         for i in model.all_species])
        return dndt

    def ic():
        res = np.empty(model.n_heights, dtype='object')
        import time
        startt = time.time()

        for h in range(model.n_heights):
            n = np.array([c.density[h, 0] for c in model.all_species])
            res[h] = solve_ivp(fun, (ts[0], te[-1]), n, method='BDF', vectorized=False, args=[h],
                               t_eval=ts_, max_step=0.44, atol = 1e-3)

            import sys
            sys.stdout.write('\r' + (' ' * 22))
            sys.stdout.write("\r Height {0}. Mean time {1}".format(h, (time.time() - startt) / (h + 1)))
            sys.stdout.flush()

            if res[h].status != 0:
                print(h)
                print(res[h])
                for c in model.all_species:
                    plt.figure()
                    plt.plot(res[h].t, res[h].y[c.c_ID, :], label='n(' + c.name + ')')
                    plt.yscale('log')
                    plt.tight_layout()
                    plt.legend(loc=2)
                    plt.show()
                breakpoint()


        #        res[h] = solve_ivp(fun, (ts[0], te[-1]), n, method='BDF',vectorized=False, args = [h],
        #                           t_eval = np.arange(0, te[-1], 0.01), max_step = 0.0444)
        # for j, c in enumerate(model.all_species):
        #    c.density[h] = res[h].y[j, -1]           #replaced by below code:

        new_densities = np.array([alt.y[:, -1] for alt in res]).T
        # return new_densities
        return res

    res = ic()

    # check charge neutrality!!
    # [e,O,Op,O2,O2p,N,Np,N2,N2p,NO,NOp,H,Hp] = np.array([r.y for r in res]).swapaxes(0, 1)

    for h, i in enumerate(res):
        for c in model.all_species:
            plt.figure()
            plt.plot(i.t, i.y[c.c_ID, :], label=c.name)
            if c == model.e: plt.plot(ts_, ne[h, :], label='ElSpec ne')
            if c == model.N2: plt.plot(ts_, nN2[h], label='ElSpec N2')
            if c == model.O2: plt.plot(ts_, nO2[h], label='ElSpec O2')
            if c == model.O: plt.plot(ts_, nO[h], label='ElSpec O ')
            if c == model.NOp: plt.plot(ts_, nNOp[h], label='ElSpec NOp')
            if c == model.O2p: plt.plot(ts_, nO2p[h], label='ElSpec O2p')
            # if c == model.Op: plt.plot(ts_, nOp[h], label='ElSpec Op ')
            plt.legend(loc=2)
            plt.yscale('log')
            plt.xlabel('Time [s]')
            plt.ylabel('Density [m-3]')
            plt.title(c.name + ' Density')
            ax2 = plt.gca().twinx()
            ax2.plot(ts_, e_prod[h, :], '.', color='green', label='q_e')
            ax2.set_yscale('log')
            ax2.legend(loc=1)
            ax2.set_ylabel('Electron production [m-3s-1]')
            # for t in ts_: plt.axvline(t, alpha = 0.1)

            plt.savefig(direc + 'plots/' + c.name + ' Density' + '_IC_' + str(iteration) + '.svg')
        break

    n_ic_ = np.array([r.y for r in res])
    n_ic = np.empty(n_ic_.shape)

    print('mixf =', mixf)
    if iteration == 0:
        if n_ic_.shape[2] == ne.shape[1]:
            print('untested')
            # breakpoint()
            n_ic_old = np.array([c.density for c in model.all_species]).swapaxes(0, 1)
            n_ic = (n_ic_ + mixf * n_ic_old) / (1 + mixf)
        else:
            # interpolate
            n_ic_old = np.array(
                [[PchipInterpolator(ts, d)(ts_int) for d in c.density] for c in model.all_species]).swapaxes(0, 1)
            n_ic = (n_ic_ + mixf * n_ic_old) / (1 + mixf)
    else:
        with open(direc + "IC_res_" + str(iteration - 1) + '.pickle', 'rb') as pf:
            [ts_, z, n_ic_old, eff_rr_old] = pickle.load(pf)
        n_ic = (n_ic_ + mixf * n_ic_old) / (1 + mixf)

    c_order = np.array([c.name for c in model.all_species])
    order = ','.join(c_order).replace('+', 'p').replace('-', '').replace('(', '_').replace(')', '')
    for c, n in zip(order.split(','), n_ic.swapaxes(0, 1)):
        #print(c)
        exec(f"global {c}; {c} = n")  # , {"n": n, f"{c}": c})
    # print(Op_4S)

    eff_rr = (rrate.T[:, 0, :] * n_ic[:, 3, :] + \
              rrate.T[:, 1, :] * n_ic[:, 6, :] + \
              rrate.T[:, 2, :] * n_ic[:, 8, :] + \
              rrate.T[:, 3, :] * n_ic[:, 15, :]) / n_ic[:, 0, :]

    [Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
    # [e, O, Op, O2, O2p, N, Np, N2, N2p, NO, NOp, H, Hp] = n_ic.swapaxes(0, 1)
    elspec_iri_sorted = np.array([Tn, Ti, Te, nN2, nO2, nO, nAr, NOp, O2p, Op_4S]).swapaxes(0, 1)

    #if iteration == 0:
    #    elspec_iri_sorted = np.array([Tn, Ti, Te, nN2, nO2, nO, nAr, NOp, O2p, Op_4S]).swapaxes(0, 1)

    mdict = {"elspec_iri_sorted": elspec_iri_sorted, "eff_rr": eff_rr}
    spio.savemat(direc + 'IC_' + str(iteration) + '.mat', mdict)

    savedir = direc + "IC_res_" + str(iteration) + '.pickle'
    print(savedir)
    with open(savedir, "wb") as f:
        pickle.dump([ts_, z_model, n_ic, eff_rr], f, protocol=pickle.HIGHEST_PROTOCOL)

    d_effrr = con["alpha"] - eff_rr

    plt.figure()
    plt.plot(ts_, np.abs(d_effrr[0, :]), label='d_alpha')
    plt.plot(ts_, eff_rr[0, :], label='alpha')
    plt.xlabel('Time [s]')
    plt.ylabel('Eff. Recombination Rate [m3s-1]')
    plt.yscale('log')
    plt.legend()
    plt.title('Eff.Rec.Rate and Difference from last Iteratrion')
    plt.savefig(direc + 'plots/IC_' + str(iteration) + '_' + \
                'Eff Rec Rate and Difference from last Iteratrion.svg')

    h = 20
    [re, rO, rOp, rO2, rO2p, rN, rNp, rN2, rN2p, rNO, rNOp, rH, rHp] = [e, O, Op_4S, O2, O2p, N_4S, Np, N2, N2p, NO,
                                                                        NOp,
                                                                        H, Hp] / e
    plt.figure()
    plt.stackplot(ts_, rOp[h], rO2p[h], rNp[h], rN2p[h], rNOp[h], rHp[h], \
                  labels=['O+', 'O2+', 'N+', 'N2+', 'NO+', 'H+'])
    plt.xlabel('Time [s]')
    plt.ylabel('Ratio of Charged Species')
    plt.legend(loc=2)
    plt.title('Charged Species Stackplot at height index ' + str(h))
    plt.savefig(direc + 'plots/IC_' + str(iteration) + '_' + \
                'Charged Species Stackplot at height index ' + str(h) + '.svg')

    # check charge nuetrality!!
    sum_charged = np.sum(np.array([Op_4S, Op_2D, Op_2P, O2p_a4P, O2p, Np, N2p, NOp, Hp]), axis=0)
    r = np.abs(sum_charged - e)
    r.shape
    plt.figure()
    pm = plt.pcolor(ts_, z_model, r / e)
    plt.gcf().colorbar(pm)
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [km]')
    plt.title('Relative Charge imabalance $(n_{I^+} - n_{e^-})/n_{e^-}$')
    plt.savefig(direc + 'plots/IC_' + str(iteration) + '_' + \
                'Relative Charge imabalance.svg')

    plt.figure()
    pc = plt.pcolor(ts_, z_model, eff_rr, label='alpha')
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [km]')
    plt.legend()
    # plt.tight_layout()
    plt.colorbar(pc)
    plt.title('Effective Recombination rate [m3s-1]')
    plt.savefig(direc + 'plots/IC_' + str(iteration) + '_' + \
                'Effective Recombination rate.svg')

    import matplotlib as mpl
    plt.figure()
    pc = plt.pcolor(ts_, z_model, d_effrr, norm=mpl.colors.CenteredNorm(),
                    label='alpha', cmap='RdBu')
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [km]')
    # plt.tight_layout()
    plt.colorbar(pc)
    plt.title('Deviation from previous eff. rec. rate [m3s-1]')
    plt.savefig(direc + 'plots/IC_' + str(iteration) + '_' + \
                'Deviation from previous eff rec rate.svg')

