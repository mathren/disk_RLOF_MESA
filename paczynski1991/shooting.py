import numpy as np
from scipy.optimize import bisect, root
from scipy.integrate import solve_ivp

# Choose units
# G = gravitational constant = 1
# M = accretor mass = 1
# Rp = accretor polar radius = 1 (this defines Rp = R_star in P+91)

# Dimensionless free parameters (Paczynski+1991)
a_inv = 1
ω_s = 0.5

n = 1.5 # polytropic index (n=3/2 is γ=5/3)

# Chosen for numerical reasons (x -> non-dim'd r)
x_out = 100 # outer boundary condition (where disk is approximately Keplerian)
x_mid = 50 # matching point
x_in = 0.8 # inner boundary condition (within star but not totally)

# Form a vector Z = [ω, y] and do shooting to find j_star
# ω = dimensionless disk rotation / orbital rate
# y = (thin) disk scale height

# Define a function giving dZ / dx:
def get_dZ_dx(Z, x, j_star):
    ω, y = Z
    dω_dx = x ** (3. * n - 1.5) * (j_star - x ** 2 * ω) * a_inv / y ** (2. * n + 3.)
    dy_dx = (x ** 2 + y ** 2) ** 1.5 * x * ω ** 2 / y - 1
    dZ_dx = np.array([dω_dx, dy_dx])

# Define a function giving boundary condition at inner Z
def get_Z_in():
    ω_in = ω_s
    y_in = ((1 - 0.5 * ω_s ** 2 * x_in ** 2) ** -2. - x_in ** 2) ** 0.5 # P91, Eqn. 27
    Z_in = np.array([ω_in, y_in])

    return Z_in

def get_Z_out(j_star):
    # Fix ω and dω/dx to be Keplerian
    ω_out = x_out ** -1.5
    y_out = (2. * (x_out ** 0.5 - j_star) * a_inv * x_out ** 2.5 / 3.) ** (1 / (2. * n + 3.))
    Z_out = np.array([ω_out, y_out])

    return Z_out

def get_Z_match_from_in(j_star):
    """
    Integrate from inner boundary condition to get Z_match
    """
    ...

def get_Z_match_from_out(j_star):
    """
    Integrate from outer boundary condition to get Z_match
    """
    ...

# Also get determinant
...

# Root solve determinant
...






