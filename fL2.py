import numpy as np
import matplotlib.pyplot as plt

"""
Binary Mass Transfer: L2 mass loss fraction model

Philip Mocz (2025)

Calculate the fraction of mass lost in the outer dick based on model by Lu+22
"""

# Physical constants (cgs units)
gravitational_constant = 6.6743e-8  # [cm^3 g^-1 s^-2]
solar_mass = 1.989e33  # [g]
solar_radius = 6.96e10  # [cm]
year = 31557600.0  # [s]
a_rad = 7.5646e-15  # radiation density constant [erg cm^-3 K^-4]
k_B = 1.3807e-16  # Boltzmann constant [erg/K]
proton_mass = 1.673e-24  # [g]
speed_of_light = 2.99792e10  # [cm/s]
mu_gas = (
    1.3 / 2.4
)  # fully ionized solar abundance disk gas (with He-to-H number ratio of 0.1)
# in general: 1 / (2 * X + 3 * Y / 4 + Z / 2) for fully ionized


def kap(rho, T):
    """
    Calculate Kramer's opacity (approximate).
    Input: rho [g/cm^3], T [Kelvin]
    Output: kappa [cm^2/g]
    """
    kappa = 0.34 + 3.0e24 * rho * T**-3.5
    return kappa


def calculate_L2_mass_loss_fraction(
    donor_mass, accretor_mass, mass_transfer_rate, orbital_separation, alpha
):
    """Calculate the L2 mass-loss fraction."""

    eps_small = 1e-12  # a very small number
    tol = 1e-8  # fractional tolerance for bisection method

    # key parameters
    M2 = accretor_mass * solar_mass
    M1dot = mass_transfer_rate * solar_mass / year
    a = orbital_separation * solar_radius
    q = accretor_mass / donor_mass  # mass ratio M2/M1

    log_q = np.log10(q)

    # positions of Lagrangian points (based on analytic fits)
    xL1 = -0.0355 * log_q**2 + 0.251 * abs(log_q) + 0.500  # [a = SMA]
    xL2 = 0.0756 * log_q**2 - 0.424 * abs(log_q) + 1.699  # [a]

    if log_q > 0:  # m2 is more massive
        xL1 = 1 - xL1
        xL2 = 1 - xL2
    mu = q / (1 + q)
    # outer disk radius
    Rd_over_a = (1 - xL1) ** 4 / mu
    # relavent potential energies
    PhiL1_dimless = -(
        (1 - mu) / abs(xL1) + mu / abs(1 - xL1) + 0.5 * (xL1 - mu) ** 2
    )  # [G(M1+M2)/a]
    PhiL2_dimless = -(
        (1 - mu) / abs(xL2) + mu / abs(1 - xL2) + 0.5 * (xL2 - mu) ** 2
    )  # [G(M1+M2)/a]
    PhiRd_dimless = -(1 - mu + mu / Rd_over_a + 0.5 * (1 - mu) ** 2)

    GM2 = gravitational_constant * M2

    Rd = Rd_over_a * a
    Phi_units = gravitational_constant * (M2 / mu) / a
    PhiL1 = PhiL1_dimless * Phi_units
    PhiL2 = PhiL2_dimless * Phi_units
    PhiRd = PhiRd_dimless * Phi_units
    # Keplerian frequency at Rd
    omega_K = np.sqrt(GM2 / Rd**3)

    # constants involved in numerical solutions
    c1 = 2 * np.pi * a_rad * alpha * Rd / (3 * omega_K * M1dot)
    c2 = k_B * Rd / (GM2 * mu_gas * proton_mass)
    c3 = 8 * np.pi**2 * a_rad * alpha * speed_of_light * Rd**2 / (M1dot**2 * omega_K)
    c4 = (
        2
        * np.pi
        * mu_gas
        * a_rad
        * alpha
        * omega_K
        * proton_mass
        * Rd**3
        / (k_B * M1dot)
    )

    # below is for grid search for disk thickness [only used at the beginning]
    n_the = 50  # ~50 is accurate enough
    the_grid_min = 0.1
    the_grid_max = 1.0
    the_grid = np.logspace(
        np.log10(the_grid_min), np.log10(the_grid_max), n_the, endpoint=True
    )
    T_floor = 3e3  # [K] --- the minimum value for disk temperature solution

    # only T < Tmax is possible to calculate
    Tmax = (4.0 / (27.0 * c1**2 * c2)) ** (1.0 / 9.0)

    # helper function
    def f1_the_T_fL2(the, T, fL2):
        return c1 * T**4 * the**3 / (1 - fL2) - the**2 + c2 * T

    Tarr = np.zeros(n_the, dtype=float)
    for i in range(n_the):
        the = the_grid[i]
        # use bisection method
        T_left = 0.1 * min(the**2 / c2, Tmax)
        f1_left = f1_the_T_fL2(the, T_left, fL2=0)
        T_right = Tmax
        while abs((T_left - T_right) / T_right) > tol:
            T = (T_left + T_right) / 2
            f1 = f1_the_T_fL2(the, T, fL2=0)
            if f1 * f1_left > 0:
                T_left = T
                f1_left = f1
            else:
                T_right = T
        Tarr[i] = (T_left + T_right) / 2
    # now we have obtained numerical relation between the and T
    logTarr = np.log10(Tarr)
    logthe_grid = np.log10(the_grid)
    dlogthe = logthe_grid[1] - logthe_grid[0]

    # helper functions
    def f2_the_T_fL2(the, T, fL2):
        x = c4 * (T * the) ** 3 / (1 - fL2)
        U_over_P = (1.5 + x) / (1 + 1.0 / 3 * x)
        rho = (1 - fL2) * M1dot / (2 * np.pi * alpha * omega_K * Rd**3) / the**3
        return (
            7.0 / 4
            - (1.5 * U_over_P + c3 * T**4 / kap(rho, T) / (1 - fL2) ** 2) * the**2
            - PhiRd / (GM2 / Rd)
            + (PhiL1 - fL2 * PhiL2) / (GM2 / Rd) / (1 - fL2)
        )

    def T_the_nofL2(the):  # under the assumption fL2=0
        # return ((the**2/c2)**(-s) + ((c1*the)**-0.25)**(-s))**(-1./s)  # analytic (not perfect)
        logthe = np.log10(the)
        if logthe > logthe_grid[-2]:  # use analytic extrapolation
            return 10 ** (logTarr[-2] - 0.25 * (logthe - logthe_grid[-2]))
        if logthe < logthe_grid[0]:  # analytic extrapolation
            return 10 ** (logTarr[0] + 2 * (logthe - logthe_grid[0]))
        i_the = int(np.floor((logthe - logthe_grid[0]) / dlogthe))
        slope = (logTarr[i_the + 1] - logTarr[i_the]) / dlogthe
        logT = logTarr[i_the] + (logthe - logthe_grid[i_the]) * slope
        return 10**logT

    def T_the(the):  # only for non-zero fL2
        return (8 * the**2 + 1 - 2 * (PhiL2 - PhiRd) / (GM2 / Rd)) / (3 * c2)

    def fL2_the(the):  # only for non-zero fL2
        T = T_the(the)
        return 1 - c1 * T**4 * the**3 / (the**2 - c2 * T)

    # bisection to find the numerical solution to f2(the, T, fL2=0)=0
    the_right = 1.0
    f2_right = f2_the_T_fL2(the_right, T_the_nofL2(the_right), fL2=0)
    separation_factor = 0.95
    the_left = separation_factor * the_right
    f2_left = f2_the_T_fL2(the_left, T_the_nofL2(the_left), fL2=0)
    while f2_left * f2_right > 0:  # need to decrease the_left
        the_right = the_left
        f2_right = f2_left
        the_left *= separation_factor
        f2_left = f2_the_T_fL2(the_left, T_the_nofL2(the_left), fL2=0)
    # now the solution is between the_left and the_right
    while abs((the_left - the_right) / the_right) > tol:
        the = (the_left + the_right) / 2
        f2 = f2_the_T_fL2(the, T_the_nofL2(the), fL2=0)
        if f2 * f2_left > 0:
            the_left = the
            f2_left = f2
        else:
            the_right = the
    # solution
    the = (the_left + the_right) / 2
    T = T_the_nofL2(the)
    the_max = np.sqrt(
        3.0 / 8 * c2 * T + 1.0 / 4 * (PhiL2 - PhiRd) / (GM2 / Rd) - 1.0 / 8
    )

    if the < the_max:  # this is the correct solution
        fL2 = eps_small
    else:
        the_min = (
            1.0 / 2 * np.sqrt((PhiL2 - PhiRd) / (GM2 / Rd) - 1.0 / 2)
        )  # corresponding to fL2=1, T=0
        # need to find the maximum corresponding to fL2=0
        # this is given by the intersection between T_the(the), T_the_nofL2(the)
        the_left = the_min
        the_right = 1.0
        f_left = T_the(the_left) - T_the_nofL2(the_left)
        while abs((the_left - the_right) / the_right) > tol:
            the = (the_left + the_right) / 2
            f = T_the(the) - T_the_nofL2(the)
            if f * f_left > 0:
                the_left = the
                f_left = f
            else:
                the_right = the
        the_max = (the_left + the_right) / 2  # this corresponds to fL2=0

        # --- numerical solution for f2(the, T, fL2)=0 under non-zero fL2

        # -- do not use exactly the_min (corresponding to T = 0, bc. kap table breaks down)
        # -- define another the_min based on T_floor (kap table won't be a problem)
        the_min = np.sqrt(
            3.0 / 8 * c2 * T_floor + 1.0 / 4 * (PhiL2 - PhiRd) / (GM2 / Rd) - 1.0 / 8
        )
        the_left = the_min
        f2_left = f2_the_T_fL2(the_left, T_the(the_left), fL2_the(the_left))
        the_right = the_max / (1 + eps_small)
        # bisection again
        while abs((the_left - the_right) / the_right) > tol:
            the = (the_left + the_right) / 2
            f2 = f2_the_T_fL2(the, T_the(the), fL2_the(the))
            if f2 * f2_left > 0:
                the_left = the
                f2_left = f2
            else:
                the_right = the
        # solution
        the = (the_left + the_right) / 2
        fL2 = fL2_the(the)

    return fL2


# Main execution block
if __name__ == "__main__":
    # Input
    donor_mass = 20.0  # in solar masses
    accretor_mass = 10.0  # in solar masses
    log_a_range = [0.15, 3.25]  # in solar radius
    log_M1dot_range = [-5.3, -1.3]  # in solar mass / year
    alpha = 0.1  # alpha disk viscosity parameter (dimensionless)
    # X = 0.7  # hydrogen mass fraction
    # Z = 0.02  # metallicity

    # Calculate grid of models
    a_vals = np.logspace(log_a_range[0], log_a_range[1], 100)
    M1dot_vals = np.logspace(log_M1dot_range[0], log_M1dot_range[1], 100)

    L2_frac_grid = np.zeros((len(M1dot_vals), len(a_vals)))

    for i, mass_transfer_rate in enumerate(M1dot_vals):
        for j, orbital_separation in enumerate(a_vals):
            L2_frac_grid[j, i] = calculate_L2_mass_loss_fraction(
                donor_mass,
                accretor_mass,
                mass_transfer_rate,
                orbital_separation,
                alpha,
            )

    # Reproduce Fig. 2 of Lu+22
    plt.figure(figsize=(8, 6))
    X, Y = np.meshgrid(np.log10(M1dot_vals), np.log10(a_vals))
    c = plt.contourf(X, Y, L2_frac_grid, levels=50, cmap="RdBu")
    plt.colorbar(c, label="L2 Mass Loss Fraction")
    plt.xlabel("log10(Mass Transfer Rate / $M_\\odot$ yr$^{-1}$)")
    plt.ylabel("log10(Orbital Separation / $R_\\odot$)")
    plt.title("L2 Mass Loss Fraction")
    plt.clim(0, 1)
    plt.tight_layout()
    plt.show()
