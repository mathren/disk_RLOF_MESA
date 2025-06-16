from MESAreader import getSrcCol, secyer
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Plot_functions import get_iRLOF

def plot_omega_div_omega_crit(hfile, ax, **kwargs):
    src, col = getSrcCol(hfile, False, False)
    omega_ratio = src[:, col.index("surf_avg_omega_div_omega_crit")]
    n= src[:, col.index("model_number")]
    ax.axhline(1,0,1, ls='--', lw=1, zorder=0)
    ax.plot(n, omega_ratio, **kwargs)


if __name__ == "__main__":
    root = '/home/mrenzo/Runs/disk_RLOF_MESA/MESA_Setup'
    path1 = root + '/working_starting_setup/'
    path2 = root + '/working_paczynski91_noL2'

    hfile_acc1 = path1 + '/LOGS2/history.data'
    hfile_acc2 = path2 + '/LOGS2/history.data'

    fig = plt.figure()
    gs = gridspec.GridSpec(150, 150)
    ax = fig.add_subplot(gs[:,:])
    plot_omega_div_omega_crit(hfile_acc1, ax, label="default")
    plot_omega_div_omega_crit(hfile_acc2, ax, label="Paczynski")
    ax.legend()
    ax.set_ylabel(r"$\omega/\omega_\mathrm{crit}$")
    ax.set_xlabel(r"model number")
    plt.show()
