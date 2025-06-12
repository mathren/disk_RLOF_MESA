import numpy as np
import matplotlib.pyplot as plt

# interpolation parameters
x = np.linspace(0,1,1000000)
x0 = 0.4
x1 = 0.6
y1 = 0
y0 = 1
yrange = y1 - y0

steepness = 1
y = np.ones(len(x))

i_acc_disk = x < x0
i_interp = (x > x0) & (x<x1)
i_crit = x>x1

y[i_acc_disk] = y0
y[i_interp] = y0*(1-np.tanh(steepness*(x[i_interp]-x0)/(x1-x0)))
y[i_crit] = y1

fig, ax = plt.subplots()
ax.scatter(x, y, c='orange', lw=3)
ax.axvline(x0,0,1, color='r')
ax.axvline(x1,0,1, color='r')
ax.axhline(y0, 0,1)
ax.axhline(y1, 0,1)
plt.show()
