import scipy.io as spio
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
import ionChem
import numpy as np
import matplotlib.pyplot as plt
import loadMSIS
import loadmat
import pickle
import mat73

def ic(direc, chemistry_config, file, iteration, mixf = 0, test = False):
    # load content of last Elspec iteration
    f = direc + file + str(iteration)
    try:
        con = loadmat.loadmat(f)["ElSpecOut"]
    except:
        print("Exception cought")
        con = mat73.loadmat(f + ".mat")["ElSpecOut"]
    ne = con["ne"].astype('float64')
    assert ne.dtype == 'float64'
    n_model = con["iri"]
    assert n_model.dtype == 'float64'
    [Tn_, Ti_, Te_, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
    [ne_, Ti, Te, _] = con["par"].swapaxes(0, 1)
    Tn = Tn_
    Tr = (Ti + Tn) / 2

    # setting start time [0] to 0:
    ts = con["ts"] - con["ts"][0]
    te = con["te"] - con["ts"][0]
    e_prod = con["q"]
    # print("WARNING: production deprecated")
    # e_prod[:, 20:] = np.zeros(e_prod[:, 20:].shape)
    # ts[-1] = 60*60*12
    # te[-1] = 12*60*60 + 0.44

    z_model = con["h"]

    ts_ = np.copy(ts)
    #generate starting point 30 min in the past:
    ts = np.array([-30*60, *ts])
    e_prod = np.array([con["q0"], *e_prod.T]).T
    ne = np.array([con["ne0"], *ne.T]).T
    lam = lambda x: np.array([x[:, 0], *x.T]).T
    [Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = map(lam, [Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp])
    ts_show = np.arange(ts[0], te[-1], 0.01)
    ts_int = ts  #set which time array to use for integration etc.

    # from scipy.interpolate import RegularGridInterpolator
    # from matplotlib import cm
    # X, Y = np.meshgrid(ts_, z_model)
    # Z = Te[:, 1:]
    # interp = RegularGridInterpolator((ts_, z_model), Z.T, method = 'pchip')
    # X_, Y_ = np.meshgrid(np.arange(ts_[0], ts_[-1], 0.1), np.arange(z_model[0], z_model[-1], 1))
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    #
    # norm = plt.Normalize(Z.min(), Z.max())
    # colors = cm.viridis(norm(Z))
    # rcount, ccount, _ = colors.shape
    #
    # surf = ax.plot_surface(X, Y, Z, rcount=rcount, ccount=ccount,
    #                        facecolors=colors, shade=False)
    # surf.set_facecolor((0, 0, 0, 0))
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # Z_ = interp((X_, Y_))
    # norm = plt.Normalize(Z_.min(), Z_.max())
    # colors = cm.viridis(norm(Z_))
    # rcount, ccount, _ = colors.shape
    # surf = ax.plot_surface(X_, Y_, Z_, rcount=rcount, ccount=ccount,
    #                        facecolors=colors, shade=False)
    # surf.set_facecolor((0, 0, 0, 0))
    #
    # plt.show()


    # normalise charged species to fit electron density
    # necessary for every timestep, as ne changes in ElSpec, but n(ions) does not
    [nNOp, nO2p, nOp] = np.array([nNOp, nO2p, nOp]) / np.sum(np.array([nNOp, nO2p, nOp]), axis=0) * ne

    model = ionChem.ionChem(chemistry_config, z_model)

    for c in model.all_species:
        c.density = nN2 * 0
        c.prod = e_prod * 0


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
    # print('Warnign: random densities in ionChem, line 106')
    # model.Np.density   = ne * 0
    # model.Hp.density   = ne * 0
    # model.N2p.density  = ne * 0
    # model.O2p.density  = ne * 0.1
    # model.NOp.density  = ne * 0.3
    # model.Op_4S.density  = ne * 0.6

    # model.N_2D.density = 0
    # model.N_4S.density = nO / nO[0] * 1e11
    # model.NO.density   = nN2 / nN2[0] *5e11
    # model.O_1D.density = nO / nO[0] * 1e12
    # model.O_1S.density = nO / nO[0] * 1e11
    # model.Op_4S.density= 0

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

    #model.check_chargeNeutrality()

    Op_prod = e_prod * 0.56 * model.O.density / \
               (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)
    O2p_prod = e_prod * 1.00 * model.O2.density / \
               (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)
    N2p_prod = e_prod * 0.92 * model.N2.density / \
               (0.92 * model.N2.density + model.O2.density + 0.56 * model.O.density)

    # create smooth function for production
    t = np.arange(-30*60, te[-1], 0.01)
    
    
    
    #-----------------start here with julia implementation---------------------------------------#
    
    
    

    stepf_e = np.array([stepped_prod_t(e_prod, i, ts, te) for i in t])
    e_prod_smooth = PchipInterpolator(t, stepf_e)

    stepf_Op = np.array([stepped_prod_t(Op_prod, i, ts, te) for i in t])
    Op_prod_smooth = PchipInterpolator(t, stepf_Op)

    stepf_O2p = np.array([stepped_prod_t(O2p_prod, i, ts, te) for i in t])
    O2p_prod_smooth = PchipInterpolator(t, stepf_O2p)

    stepf_N2p = np.array([stepped_prod_t(N2p_prod, i, ts, te) for i in t])
    N2p_prod_smooth = PchipInterpolator(t, stepf_N2p)

    # plt.figure()
    # plt.plot(ts_, e_prod[0], 'x')
    # plt.plot(ts_show, e_prod_smooth(ts_show)[:, 0])
    # plt.yscale('log')
    # plt.ylabel('Electron Production')
    # plt.xlabel('Time')
    # plt.show()

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
        if c.name == 'O+(4S)':prodMat[i]=Op_prod_smooth
        if c.name == 'O2+': prodMat[i] = O2p_prod_smooth
        if c.name == 'N2+': prodMat[i] = N2p_prod_smooth

    # new way to do ionospheric chemistry
    # 1. take all info from reaction (reaction rate, constituents)
    ode_mat = np.zeros((len(model.all_reactions), len(model.all_species)), dtype='object')
    rrate = np.array([np.array([r.r_rate_t(*temp) for r in model.all_reactions]) for temp in zip(Tn.T, Ti.T, Te.T)])
#    rrate_smooth = np.array([PchipInterpolator(t, np.array([stepped_prod_t(r.T, i, ts, te) for i in t])) for r in rrate.swapaxes(0, 1)])

    #investigating why there is no variation in alpha:
    #we should see variation in Temp
    plt.figure()
    plt.pcolormesh(ts, z_model, rrate[:, 0, :].T)
    plt.colorbar()
    # plt.figure()
    # plt.pcolormesh(np.arange(0, 60, 0.001), z_model, rrate_smooth[0](np.arange(0, 60, 0.001)).T)
    # plt.colorbar()
    #plt.show()
    # plt.figure()
    # plt.pcolormesh(ts, z_model, Te)
    # plt.show()

    # loss_mat = np.zeros([len(model.all_reactions), len(model.all_species), len(z_model), len(ts)])
    # for i, r in enumerate(model.all_reactions):
    #     ed0, ed1 = r.educts_ID
    #     loss_mat[i, ed0, :, :] = rrate[:, i, :].T * model.all_species[ed1].density
    #     loss_mat[i, ed1, :, :] = rrate[:, i, :].T * model.all_species[ed0].density
    # lossRate = np.sum(loss_mat, axis = 0)
    # plt.figure()
    # for i, c in enumerate(model.all_species):
    #     if c.name not in ['O2+', 'NO+', 'N2+']: continue
    #     line = plt.plot(lossRate[i, :, 0], z_model / 1e3, label=c.name)
    #     plt.text(lossRate[i, 0, 0], z_model[1] / 1e3, c.name, color=line[0].get_color())
    # plt.xscale('log')
    # plt.xlabel('Loss rate [m-3 s-1]')
    # plt.ylabel('Altitude [km]')
    # plt.legend()
    #plt.show()


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
        if i.name in ['H', 'O', 'O2', 'N2']:
            ode_raw[i.c_ID] = np.zeros([1, 4])
    #print(ode_raw)

    def asint(arr):
        return arr.astype(int)

    # 3. produce DG with only relevant terms
    def fun(t, n, h, temp):
        rrate = np.array([r.r_rate_t2(*temp(t)) for r in model.all_reactions])

        dndt = np.array([np.sum(((rrate[asint(ode_raw[i.c_ID][:, 0])].T * ode_raw[i.c_ID][:, 1]).T \
                                 * n[asint(ode_raw[i.c_ID][:, 2])] * n[asint(ode_raw[i.c_ID][:, 3])] \
                                 ), axis=0) \
                         + prodMat[i.c_ID](t)[h] \
                         for i in model.all_species])
        return dndt

    def solve_ic():
        res = np.empty(model.n_heights, dtype='object')
        import time
        startt = time.time()

        for h in range(model.n_heights):
            #breakpoint()
            n = np.array([c.density[h, 0] for c in model.all_species])
            #print([c.name for c in model.all_species])
            t = np.arange(-30*60, te[-1], 0.01)
            temp = PchipInterpolator(ts, np.array([Tn[h, :], Ti[h, :], Te[h, :]]).T)

            res[h] = solve_ivp(fun, (ts[0], te[-1]), n, method='BDF', vectorized=False, args=[h, temp],
                               t_eval=ts_int, max_step=0.44, atol = 1e-3, rtol = 1e-7, dense_output=True)

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

    res = solve_ic()

    # check charge neutrality!!
    # [e,O,Op,O2,O2p,N,Np,N2,N2p,NO,NOp,H,Hp] = np.array([r.y for r in res]).swapaxes(0, 1)

    n_ic_ = np.array([r.y for r in res])

    if any(n_ic_.flat < 0):
        print('Negative ion densities detected.')
        for i, c in enumerate(n_ic_.swapaxes(1, 2).copy()):
            if any(c.flat < 0):
                print('Negative ion densities in ', i)

    print('mixf =', mixf)
    if mixf == 0:
        n_ic = n_ic_
    else:
        if iteration == 0:
            if n_ic_.shape[2] == ne.shape[1]:
                print('untested')
                # breakpoint()
                n_ic_old = np.array([c.density for c in model.all_species]).swapaxes(0, 1)
            else:
                #interpolate
                n_ic_old = np.array([[PchipInterpolator(ts, d)(ts_int) for d in c.density] for c in model.all_species]).swapaxes(0, 1)
        else:
            with open(direc + "IC_res_" + str(iteration - 1) + '.pickle', 'rb') as pf:
                data = pickle.load(pf)
                n_ic_old = data[2]
        n_ic = (n_ic_ + mixf * n_ic_old) / (1 + mixf)

    c_order = np.array([c.name for c in model.all_species])
    order = ','.join(c_order).replace('+', 'p').replace('-', '').replace('(', '_').replace(')', '')
    for c, n in zip(order.split(','), n_ic.swapaxes(0, 1)):
        print(c)
        exec(f"global {c}; {c} = n")  # , {"n": n, f"{c}": c})

    # print(Op_4S)
    ne_init = e[:, 1]

    eff_rr = (rrate.T[:, 0, :] * n_ic[:, 3, :] + \
              rrate.T[:, 1, :] * n_ic[:, 6, :] + \
              rrate.T[:, 2, :] * n_ic[:, 8, :] + \
              rrate.T[:, 3, :] * n_ic[:, 15, :]) / n_ic[:, 0, :]

    #[Tn, Ti, Te, nN2, nO2, nO, nAr, nNOp, nO2p, nOp] = n_model.swapaxes(0, 1)
    # [e, O, Op, O2, O2p, N, Np, N2, N2p, NO, NOp, H, Hp] = n_ic.swapaxes(0, 1)
    elspec_iri_sorted = np.array([Tn_, Ti_, Te_, nN2[:, 1:], nO2[:, 1:], nO[:, 1:], nAr[:, 1:], NOp[:, 1:], O2p[:, 1:], Op_4S[:, 1:]]).swapaxes(0, 1)

    # if iteration == 0:
    #    elspec_iri_sorted = np.array([Tn, Ti, Te, nN2, nO2, nO, nAr, NOp, O2p, Op_4S]).swapaxes(0, 1)

    mdict = {"elspec_iri_sorted": elspec_iri_sorted, "eff_rr": eff_rr[:, 1:], "ne_init": ne_init}
    spio.savemat(direc + 'IC_' + str(iteration) + '.mat', mdict)

    return 0

    savedir = direc + "IC_res_" + str(iteration) + '.pickle'
    print(savedir)
    with open(savedir, "wb") as f:
        pickle.dump([ts_int, z_model, n_ic, eff_rr, ne_init, res], f, protocol=pickle.HIGHEST_PROTOCOL)

    #breakpoint()
    for h, i in enumerate(res):
        for c in model.all_species:
            plt.figure()
            ax = plt.gca()
            plt.plot(ts_show, i.sol(ts_show)[c.c_ID, :], label=r'IC n('+c.name+')')
            #plt.plot(i.t, i.y[c.c_ID, :], label=r'$n_{'+c.name+'}$')
            if c == model.e:
                plt.plot(ts_int, ne[h, :], label=r'ElSpec $n_e$')
                ax.plot(ts_int, np.sqrt(e_prod_smooth(ts_int)[:, 0] / eff_rr[0, :]), '--', color='red', label=r'$ne_{ss}$')
                ax.plot(ts_show, i.sol(ts_show)[c.c_ID, :], label = r'$ne_{ode}$')
            if c == model.N2: plt.plot(ts_, nN2[h, 1:], label='ElSpec n(N2)')
            if c == model.O2: plt.plot(ts_, nO2[h, 1:], label='ElSpec n(O2)')
            if c == model.O: plt.plot(ts_ , nO[h, 1:],  label='ElSpec n(O) ')
            if c == model.NOp: plt.plot(ts_, nNOp[h, 1:], label='ElSpec n(NO+)')
            if c == model.O2p: plt.plot(ts_, nO2p[h, 1:], label='ElSpec n(O2+)')
            # if c == model.Op: plt.plot(ts_, nOp[h], label='ElSpec n(O+) ')
            #ax.plot(ts_int, np.sqrt(e_prod_smooth(ts_int)[:, 0]/eff_rr[0, :]), '--', color = 'red', label = 'ne_ss')
            plt.legend(loc=2)
            plt.yscale('log')
            plt.xlabel('Time [s]')
            plt.ylabel(r'Density [m$^{-3}$]')
            #plt.title(c.name + ' Density')
            ax = plt.gca()
            #ax3 = ax.twinx()
            #temp = PchipInterpolator(ts, np.array([Tn[h, :], Ti[h, :], Te[h, :]]).T)
            #ax3.plot(ts_show, temp(ts_show)[:, 2], "--", label = "Te")
            #ax3.legend()
            ax2 = ax.twinx()
            ax2.plot(ts, e_prod[h, :], '.', color='green', label=r'$q_e$')
            if c == model.e:
                ax2.plot(ts_show, e_prod_smooth(ts_show)[:, 0], '--', color = 'green', label = r'$q_{e, intp}$')
            ax2.set_yscale('log')
            ax2.legend(loc=1)
            ax2.set_ylabel(r'Electron production [m$^{-3}$s$^{-1}$]')
            # for t in ts_: plt.axvline(t, alpha = 0.1)
            plt.tight_layout()
            if test == True:
                plt.show()
            plt.savefig(direc + 'plots/' + c.name + '_Density_IC_' + str(iteration) + '.svg')
            plt.savefig(direc + 'plots/' + c.name + '_Density_IC_' + str(iteration) + '.eps')
        break


    d_effrr = con["alpha"] - eff_rr[:, 1:]

    plt.figure()
    plt.plot(ts_, np.abs(d_effrr[0, :]), label='d_alpha')
    plt.plot(ts_, eff_rr[0, 1:], label='alpha')
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
    plt.stackplot(np.array([-1, *ts_]), rOp[h], rO2p[h], rNp[h], rN2p[h], rNOp[h], rHp[h], \
                  labels=['O+', 'O2+', 'N+', 'N2+', 'NO+', 'H+'])
    plt.xlabel('Time [s]')
    plt.ylabel('Ratio of Charged Species')
    plt.legend(loc=2)
    plt.title('Charged Species Stackplot at height index ' + str(h))
    plt.savefig(direc + 'plots/Charged Species Stackplot at height index ' + str(h) + 'IC_' + str(iteration) + '.svg')

    # check charge nuetrality!!
    sum_charged = np.sum(np.array([Op_4S, Op_2D, Op_2P, O2p_a4P, O2p, Np, N2p, NOp, Hp]), axis=0)
    r = np.abs(sum_charged - e)
    r.shape
    plt.figure()
    pm = plt.pcolor(np.array([-1, *ts_]), z_model, r / e)
    plt.gcf().colorbar(pm)
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [km]')
    plt.title('Relative Charge imabalance $(n_{I^+} - n_{e^-})/n_{e^-}$')
    plt.savefig(direc + 'plots/Relative Charge imbalance IC_' + str(iteration) +'.svg')

    plt.figure()
    pc = plt.pcolor(np.array([-1, *ts_]), z_model, eff_rr, label='alpha')
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [km]')
    plt.legend()
    # plt.tight_layout()
    plt.colorbar(pc)
    plt.title('Effective Recombination rate [m3s-1]')
    plt.savefig(direc + 'plots/Effective Recombination rate IC_' + str(iteration)+ '.svg')

    import matplotlib as mpl
    plt.figure()
    pc = plt.pcolor(ts_, z_model, d_effrr, norm=mpl.colors.CenteredNorm(),
                    label='alpha', cmap='RdBu')
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [km]')
    # plt.tight_layout()
    plt.colorbar(pc)
    plt.title('Deviation from previous eff. rec. rate [m3s-1]')
    plt.savefig(direc + 'plots/Deviation from previous eff rec rate IC_' + str(iteration) + '.svg')

    plt.close('all')

def stepped_prod_t(prod, t, ts, te):
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


def stepped_rrate_t(r, t, ts, te):
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
        prod_t = r[i_max_ts]
        return prod_t
