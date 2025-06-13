import numpy as np
import matplotlib.pyplot as plt

# interpolation parameters
x = np.linspace(0,1,1000)
x0 = 0.75
x1 = 0.9
y0 = 1e31
y1 = 0.0
steepness = 10

i_acc_disk = x < x0
i_interp = (x >= x0) & (x<=x1)
i_crit = x>x1

y = np.ones(len(x))
y[i_acc_disk] = y0
y[i_interp] = y0+0.5*(y1-y0)*(1+np.tanh(steepness*(x[i_interp]-(x1+x0)/2.0)/(x1-x0)))
y[i_crit] = y1

fig, ax = plt.subplots()

ax.scatter(x, y, c='orange', lw=3)
ax.axvline(x0,0,1, color='g')
ax.axvline(x1,0,1, color='r')
ax.axhline(y0, 0,1)
ax.axhline(y1, 0,1)
ax.set_yscale('log')
# ax.plot(x, (1-np.tanh(steepness*x))/2.0)
plt.show()
