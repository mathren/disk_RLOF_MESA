from MESAreader import getSrcCol, secyer
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Plot_functions import get_iRLOF



def plot_sep(hfile, ax, bx=None, **kwargs):
    iRLOF = get_iRLOF(hfile)

    src, col = getSrcCol(hfile)
    sep = src[iRLOF, col.index("binary_separation")]
    t = src[iRLOF, col.index("age")]/1e6   # Myr
    ax.plot(t, sep, **kwargs)


def plot_P(hfile, ax, bx=None, **kwargs):
    iRLOF = get_iRLOF(hfile)

    src, col = getSrcCol(hfile)
    per = src[iRLOF, col.index("period_days")]
    t = src[iRLOF, col.index("age")]/1e6   # Myr
    ax.plot(t, per, **kwargs)


def plot_masses(hfile, ax, **kwargs):
    iRLOF = get_iRLOF(hfile)

    src, col = getSrcCol(hfile)
    m1 = src[iRLOF, col.index("star_1_mass")]
    m2 = src[iRLOF, col.index("star_2_mass")]
    t = src[iRLOF, col.index("age")]/1e6   # Myr
    ax.plot(t, m1, **kwargs)
    ax.plot(t, m2, **kwargs)


def plot_checks(hfile, bhfile, ax, **kwargs):
    """ plot things between 0 and 1 """
    iRLOF = get_iRLOF(bhfile)
    src, col = getSrcCol(bhfile)
    t = src[iRLOF, col.index("age")]/1e6   # Myr
    m1 = src[iRLOF, col.index("star_1_mass")]
    m2 = src[iRLOF, col.index("star_2_mass")]
    accretion_mode = src[iRLOF, col.index("accretion_mode")]
    try:
        src, col = getSrcCol(hfile)
        surf_avg_omega_div_omega_crit = src[iRLOF, col.index("surf_avg_omega_div_omega_crit")]
        ax.plot(t, surf_avg_omega_div_omega_crit, c='C1', **kwargs)
    except IndexError:
        pass
    ax.plot(t, (m1+m2)/(m1[0]+m2[0]), c='C0', **kwargs)
    ax.plot(t, accretion_mode, c='C2', **kwargs)

if __name__ == "__main__":
    root = '/home/mrenzo/Runs/disk_RLOF_MESA/MESA_Setup'
    path1 = root + '/working_starting_setup/'
    path2 = root + '/working_paczynski91_noL2'
    path3 = root
    binhfile1 = path1 + "/binary_history.data"
    binhfile2 = path2 + "/binary_history.data"
    binhfile3 = path3 + "/binary_history.data"

    hfile_acc1 = path1 + '/LOGS2/history.data'
    hfile_acc2 = path2 + '/LOGS2/history.data'
    hfile_acc3 = path3 + '/LOGS2/history.data'

    label1 = "no L2, yes jdot"
    label2 = "no L2, Paczynski91"
    label3 = "L2, Paczynski91"

    fig = plt.figure()
    gs = gridspec.GridSpec(200, 150)

    ax = fig.add_subplot(gs[:50,:])
    bx = fig.add_subplot(gs[50:100,:])
    cx = fig.add_subplot(gs[100:150,:])
    dx = fig.add_subplot(gs[150:,:])

    dx.set_xlabel(r"age [Myr]")

    ax.set_ylabel(r"Separation [$R_\odot$]")
    bx.set_ylabel(r"P [days]")
    cx.set_ylabel(r"$M \ [M_\odot]$")
    dx.set_ylabel(r"checks")

    plot_sep(binhfile1, ax, ls='-', label=label1)
    plot_sep(binhfile2, ax, ls='--', lw=4, label=label2)
    plot_sep(binhfile3, ax, ls='-.', lw=5, label=label3)
    plot_P(binhfile1, bx, ls='-', label=label1)
    plot_P(binhfile2, bx, ls='--', lw=4, label=label2)
    plot_P(binhfile3, bx, ls='-.', lw=5, label=label3)
    plot_masses(binhfile1, cx, ls='-', label=label1)
    plot_masses(binhfile2, cx, ls='--', lw=4, label=label2)
    plot_masses(binhfile3, cx, ls='-.', lw=5, label=label3)
    plot_checks(hfile_acc1, binhfile1, dx, ls='-', label=label1)
    plot_checks(hfile_acc2, binhfile2, dx, ls='--', lw=4, label=label2)
    plot_checks(hfile_acc3, binhfile3, dx, ls='-.', lw=5, label=label3)

    bx.legend()
    bx.set_xlim(ax.get_xlim())
    cx.set_xlim(ax.get_xlim())
    plt.show()
