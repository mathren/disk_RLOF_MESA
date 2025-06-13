from MESAreader import getSrcCol, secyer
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def get_total_AM(bin_hfile, single_hfile):
    src, col = getSrcCol(bin_hfile)
    src_single, col_single = getSrcCol(single_hfile)
    J_orb = src[:, col.index("J_orb")]
    J_1 = src[:, col.index("J_spin_1")]
    J_2 = src[:, col.index("J_spin_2")]
    J_lost = src[:, col.index("Jdot")]*(10**(src_single[:, col_single.index("log_dt")]))
    J_tot = J_1+J_2+J_orb+J_lost
    return J_tot


def get_total_mass(bin_hfile):
    src, col = getSrcCol(bin_hfile)
    m1 = src[:, col.index("star_1_mass")]
    m2 = src[:, col.index("star_2_mass")]
    dt = 10.0**(src[:, col.index("log_dt")]) # sec
    wind1 = 10.0**(src[:, col.index("lg_wind_mdot_1")])*dt/secyer
    wind2 = 10.0**(src[:, col.index("lg_wind_mdot_2")])*dt/secyer
    fL2 = src[:, col.index("fL2")]
    L2 = fL2*(10.0**(src[:, col.index("lg_mtransfer_rate")])*dt/secyer)
    return m1+m2+wind1+wind2+L2 # should be constant

if __name__ == "__main__":
    fig = plt.figure()
    gs = gridspec.GridSpec(150, 100)
    ax = fig.add_subplot(gs[:, :50])
    bx = ax.twinx()
    root = "../MESA_Setup/20_18_100/"
    bin_hfile = root+'binary_history.data'
    h1 = root+'LOGS1/history.data'
    h2 = root+'LOGS2/history.data'
    src, col = getSrcCol(bin_hfile)
    fL2 = src[:, col.index("fL2")]
    MdotL2 = src[:, col.index("mdot_L2")]
    accretion_mode = src[:, col.index("accretion_mode")]
    time = src[:, col.index("model_number")]# src[:, col.index("age")]
    ax.plot(time, fL2)
    ax.plot(time, mdot_L2)
    bx.plot(time, accretion_mode, c='C2')
    bx.set_ylabel(r"Accretion mode", color='C2')
    ax.set_ylabel(r"$\log_{10}(f_{L2})$")
    ax.set_yscale('log')
    # J_tot = get_total_AM(bin_hfile, h1)

    # rl_relative_overflow_1 = src[:, col.index("rl_relative_overflow_1")]
    # print(min(rl_relative_overflow_1))
    # iRLOF = (rl_relative_overflow_1 < 0)
    # ax.axvline(min(time[iRLOF]))
    # ax.axvline(max(time[iRLOF]))
    # ax.plot(time, J_tot)

    plt.show()
