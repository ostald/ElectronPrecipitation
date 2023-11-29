import matplotlib.pyplot as plt
def plot_nop_o2p_prod(N2p_prod, O2p_prod, Op_prod, model, n_ic, z_model, temp):
    ii = 283
    plt.figure()
    for r in model.all_reactions:
        if "NO+" in r.products:
            plt.plot(r.r_rate_t(*temp[:, :, ii]) * n_ic[:, r.educts_ID[0], ii] * n_ic[:, r.educts_ID[1], ii], z_model,
                     label=r.r_stoch)
    for r in model.all_reactions:
        if "O2+" in r.products:
            plt.plot(r.r_rate_t(*temp[:, :, ii]) * n_ic[:, r.educts_ID[0], ii] * n_ic[:, r.educts_ID[1], ii], z_model,
                     label=r.r_stoch, linestyle='--')
    plt.plot(Op_prod[:, ii], z_model, label='e- + O => 2e- + O+', linestyle=':')
    plt.plot(O2p_prod[:, ii], z_model, label='e- + O2 => 2e- + O2+', linestyle=':')
    plt.plot(N2p_prod[:, ii], z_model, label='e- + N2 => 2e- + N2+', linestyle=':')
    plt.xscale("log")
    plt.legend()
    plt.show()
    return 0