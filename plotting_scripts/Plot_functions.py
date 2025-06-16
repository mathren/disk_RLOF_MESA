from MESAreader import getSrcCol

def get_iRLOF(hfile, broaden=1e-3):
    src, col = getSrcCol(hfile)
    rl_relative_overflow_1 = src[:, col.index("rl_relative_overflow_1")]
    # rl_relative_overflow_2 = src[:, col.index("rl_relative_overflow_2")]
    iRLOF = rl_relative_overflow_1 >= 0-broaden
    return iRLOF
