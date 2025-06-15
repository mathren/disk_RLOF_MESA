from MESAreader import getSrcCol, secyer
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def plot_omega_div_omega_crit(hfile, ax, **kwargs):
    src, col = getSrcCol(hfile)
