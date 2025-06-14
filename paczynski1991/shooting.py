import numpy as np
from scipy.optimize import bisect, root
from scipy.integrate import solve_ivp, odeint

import matplotlib.pyplot as plt

import pdb

# Choose units
# G = gravitational constant = 1
# M = accretor mass = 1
# Rp = accretor polar radius = 1 (this defines Rp = R_star in P+91)

# Dimensionless free parameters (Paczynski+1991)
log10_a = 6
# ω_s = 0.1111 #
ω_s = np.sqrt(8 / 27)-0.04      #  0.5

n = 1.5  # polytropic index (n=3/2 is γ=5/3)

# Chosen for numerical reasons (x -> non-dim'd r)
x_out = 1.23  # outer boundary condition (where disk is approximately Keplerian)
x_in = 0.8  # inner boundary condition (within star but not totally)
method = "Radau"

# Form a vector Z = [ω, y] and do shooting to find j_star
# ω = dimensionless disk rotation / orbital rate
# y = (thin) disk scale height

# Define a function giving dZ / dx:
def get_dZ_dx(x, Z, j_star):
    ω, y = Z
    # dω_dx = x ** 3 * (j_star - ω * x ** 2) * 10 ** -log10_a / y ** 6 # n=3/2 specifically
    # dy_dx = ω ** 2 * (x ** 2 + y ** 2) ** 1.5 * x / y - x / y

    # dω_dx = x ** (3. * n - 1.5) * (j_star - ω * x ** 2) * 10 ** -log10_a / y ** (2. * n + 3.)

    dω_dx = x ** 3 * (j_star - np.exp(ω) * x ** 2) * 10 ** -log10_a * np.exp(-6 * y) * np.exp(-ω) # n=3/2 specifically
    dy_dx = (np.exp(2 * ω) * (x ** 2 + np.exp(2 * y)) ** 1.5 - 1) * x * np.exp(-2 * y)

    # dω_dx = x ** 3 * (j_star - ω) * 10 ** -log10_a / y ** 6 # n=3/2 specifically
    # dy_dx = ω ** 2 * (1 + y ** 2) ** 1.5 / y - 1 / y



    dZ_dx = np.array([dω_dx, dy_dx])

    return dZ_dx

# Define a function giving boundary condition at inner Z
def get_Z_in():
    # ω_in = ω_s
    # y_in = ((1 - 0.5 * ω_s**2 * x_in**2) ** -2.0 - x_in**2) ** 0.5  # P91, Eqn. 27
    
    ω_in = np.log(ω_s)
    y_in = np.log(((1 - 0.5 * ω_s**2 * x_in**2) ** -2.0 - x_in**2) ** 0.5)  # P91, Eqn. 27
    Z_in = np.array([ω_in, y_in])

    return Z_in





# def get_Z_out(j_star):
#     # Fix ω and dω/dx to be Keplerian
#     ω_out = x_out ** -1.5
#     y_out = (2. * 10 ** -log10_a / 3.) ** (1. / (2. * n + 3.)) * x_out ** ((6. * n + 3.) / (4. * n + 6.)) * (1 - j_star * x_out ** -0.5) ** (1. / (2. * n + 3.))
#     # y_out = (2. * 10 ** -log10_a / 3.) ** (1. / 6.) * x_out * (1 - j_star * x_out ** -0.5) ** (1. / 6.)
#     Z_out = np.array([ω_out, y_out])

#     return Z_out


def get_Z_out_from_in(j_star):
    """
    Integrate from inner boundary condition to get Z_match

    TODO: seems to break for j_star near one as dZ/dx diverges... is that physical?
    """
    x_arr = np.array([x_in, x_out])

    atol = 1e-10

    res = solve_ivp(
        get_dZ_dx,
        [x_in, x_out],
        get_Z_in(),
        t_eval=[x_out],
        args=[j_star],
        method=method,
        atol=atol
    )

    assert res.success
    Z_match = res.y[:, 0]

    # test
    x_arr = np.linspace(x_in, x_out, 10000)
    res = solve_ivp(
        get_dZ_dx,
        [x_in, x_out],
        get_Z_in(),
        t_eval=x_arr,
        args=[j_star],
        method=method,
        atol=atol
    )

    x, Z = res.t, res.y
    ω, y = res.y[0], res.y[1]

    x_a = np.linspace(0, 1, 100)
    # y_a = ((1 - 0.5 * ω_s**2 * x_a**2) ** -2.0 - x_a**2) ** 0.5
    y_a = np.log(((1 - 0.5 * ω_s**2 * x_a**2) ** -2.0 - x_a**2) ** 0.5)

    plt.close()
    plt.plot(x, ω, label="ω")
    plt.plot(x, y, label="y")
    plt.plot(x_a, y_a)
    # plt.yscale("log")
    # plt.xscale('log')
    plt.legend(loc="lower left")
    plt.show()

    return Z_match


# Also get determinant
...

# Root solve determinant
...
