#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 30 13:12:16 2025

@author: peterscherbak
"""



import numpy as np
from scipy.optimize import root, minimize_scalar
import matplotlib.pyplot as plt

# === Parameters ===
n = 1.5
alpha_fixed = 1 / 0.1  # viscosity parameter
#alpha_fixed = 1 # viscosity parameter
#alpha_fixed = .1 
#alpha_fixed = 1e6
alpha_fixed = 1e6
#
r0, r1 = 0.8, 10
#r_match = 5
N = 200  # Number of grid points

omega_over_omega_crit = 1.0 #0.5 #1.0
#omega_over_omega_crit = 0.1

#r = np.linspace(r0, r1, N)
r = np.logspace(np.log10(r0), np.log10(r1), N)





#dr = r[1] - r[0]
dr_forward = np.diff(r)
dr_backward = np.roll(dr_forward, 1)
dr_center = r[2:] - r[:-2]


omega_guess = r**-1.5
z_guess = np.linspace(0.1, 0.9, N)





# === ODE Right-hand sides ===
def f1(r, omega, z, j, alpha):
    #z = max(z, 1e-4)
    return r**(3*n - 1.5) / z**(2*n + 3) * (j - r**2 * omega) / alpha

def f2(r, omega, z, j, alpha):
    #z = max(z, 1e-4)
    return r / z * ((r**2 + z**2)**(3/2) * omega**2 - 1)

# === Boundary conditions ===
def inner_bc(j):
    
    omega0 = omega_over_omega_crit * (1.5)**(-1.5)      #1.0 #1.0 #0.5 #1.0 # 0.5
    z0 = ((1 - 0.5 * omega0**2 * r0**2)**-2. - r0**2)**0.5
    return omega0, z0

def outer_bc(j):
    omega1 = r1**-1.5
    omega1 = -1.5 * r1**(-2.5)
    term1 = (2/3 / alpha_fixed)**(1 / (2*n + 3))
    term2 = r1
    val = 1 - j * r1**(-0.5)
    if val <= 0:
        raise ValueError("Invalid outer boundary: 1 - j*r^(-1/2) <= 0")
    term3 = val**(1 / (2*n + 3))
    z1 = term1 * term2 * term3
    return omega1, z1

# === Residual function for relaxation ===
def relaxation_residual(U, j, alpha):
    omega = U[:N]
    z = U[N:]
    
    res = np.zeros_like(U)
    
    # Internal points (central difference)
    for i in range(1, N-1):
        #domega_dr = (omega[i+1] - omega[i-1]) / (2*dr)
        #dz_dr = (z[i+1] - z[i-1]) / (2*dr)
        
        if False:
            dr = r[i+1] - r[i-1]
            domega_dr = (omega[i+1] - omega[i-1]) / dr
            dz_dr = (z[i+1] - z[i-1]) / dr
    
            res[i] = domega_dr - f1(r[i], omega[i], z[i], j, alpha)
            res[N+i] = dz_dr - f2(r[i], omega[i], z[i], j, alpha)
            
        elif True:
            
            eps=1e-4
            dr_forward = r[i+1] - r[i]
            dr_backward = r[i] - r[i-1]
        
            domega_dr = (omega[i+1] - omega[i-1]) / (dr_forward + dr_backward)
            dz_dr = (z[i+1] - z[i-1]) / (dr_forward + dr_backward)
        
            d2omega_dr2 = ( (omega[i+1] - omega[i]) / dr_forward - (omega[i] - omega[i-1]) / dr_backward ) / ((dr_forward + dr_backward)/2)
            d2z_dr2 = ( (z[i+1] - z[i]) / dr_forward - (z[i] - z[i-1]) / dr_backward ) / ((dr_forward + dr_backward)/2)
        
            res[i] = domega_dr - f1(r[i], omega[i], z[i], j, alpha) + eps * d2omega_dr2
            res[N+i] = dz_dr - f2(r[i], omega[i], z[i], j, alpha) + eps * d2z_dr2

        #res[i] = domega_dr - f1(r[i], omega[i], z[i], j, alpha)
        #res[N+i] = dz_dr - f2(r[i], omega[i], z[i], j, alpha)

    # Boundary conditions
    omega0, z0 = inner_bc(j)
    omega1, z1 = outer_bc(j)
    
    domega_dr_outer = (omega[-1] - omega[-2]) / (r[-1] - r[-2])
   

    res[0] = omega[0] - omega0
    #res[N-1] = omega[-1] - omega1
    res[N-1] = domega_dr_outer - omega1

    res[N] = z[0] - z0
    #res[2*N - 1] = z[-1] - z1

    return res


def residual_norm_for_j(j):
    try:
        omega0, z0 = inner_bc(j)
        omega1, z1 = outer_bc(j)
    except ValueError:
        return 1e6  # Invalid boundary, large penalty

    #omega_guess = np.linspace(omega0, omega1, N)
    #omega_guess = r**-1.5
    #z_guess = np.linspace(z0, z1, N)
    
    
    term1 = (2/3 / alpha_fixed)**(1 / (2*n + 3))
    term2 = r
    val = 1 - j * r**(-0.5)
    bad_indices = np.where(val <= 0)[0]
    if len(bad_indices) > 0:
        do_nada=True
        #print('bad val')
        #return 1e6
        
        #print(val[bad_indices])
        
    term3 = val**(1 / (2*n + 3))
    #z_guess = term1 * term2 * term3 !!!!! had this earlier
    
    #print(z0, z1)
    
    #z_guess = np.linspace(z0, z1, N)
    #z_guess = np.linspace(0.1, 0.9, N)
    
    
    
    
    U0 = np.concatenate([omega_guess, z_guess])

    sol = root(relaxation_residual, U0, args=(j, alpha_fixed), method='lm', tol=1e-3)
    if not sol.success:
        return 1e6  # Failed to converge, large penalty

    return np.linalg.norm(relaxation_residual(sol.x, j, alpha_fixed))

# === Minimize the residual norm to find best j ===
#result = minimize_scalar(residual_norm_for_j, bounds=(-10, 10), method='bounded', options={'xatol': 1e-3})


jmax_bound = 1*r1**(1/2) # !!!!probably want to go higher!!!
#jmax_bound = 3
#result = minimize_scalar(residual_norm_for_j, bounds=(0, jmax_bound), method='bounded', options={'xatol': 1e-6})
result = minimize_scalar(residual_norm_for_j, bounds=(0, jmax_bound), method='bounded', options={'xatol': 1e-6})


#%%
#omega_over_omega_crit = 1.0
if result.success:
    best_j = result.x
    print(f"✅ Found optimal j = {best_j:.6f}")

    # Solve again with best j for plotting
    omega0, z0 = inner_bc(best_j)
    omega1, z1 = outer_bc(best_j)
    #omega_guess = np.linspace(omega0, omega1, N)
    
    
    #omega_guess = r**-1.5
    #omega_guess = r**-1.5

    #z_guess = np.linspace(z0, z1, N)
    #'''
    term1 = (2/3 / alpha_fixed)**(1 / (2*n + 3))
    term2 = r
    val = 1 - best_j * r**(-0.5)
    #if val <= 0:
         #raise ValueError("Invalid outer boundary: 1 - j*r^(-1/2) <= 0")
    term3 = val**(1 / (2*n + 3))
    #'''
    #z_guess = term1 * term2 * term3
    
    #z_guess = np.linspace(0.1, 0.9, N)
    U0 = np.concatenate([omega_guess, z_guess])

    sol = root(relaxation_residual, U0, args=(best_j, alpha_fixed), method='lm', tol=1e-3)

    if sol.success:
        omega_sol = sol.x[:N]
        z_sol = sol.x[N:]
        
        fig, ax1=plt.subplots()        


        plt.plot(r, omega_sol, label=r'$\omega(r)$', color='blue')
        plt.plot(r, z_sol, label=r'$z(r)$', color='orange')
        
        plt.plot(r, omega_guess, '--', label='Initial $\omega$', color='blue', alpha=0.5)
        plt.plot(r, z_guess, '--', label='Initial $z$', color='orange', alpha=0.5)
        
        plt.scatter(r0, omega0, color='blue')
        plt.scatter(r1, omega1, color='blue')
        
        plt.scatter(r0, z0, color='orange')
        #plt.scatter(r1, z1, color='orange')

        
        
        plt.xlabel('r')
        #plt.title(f'Solution with optimal j = {best_j:.6f} for alpha = {alpha_fixed:.2f}, $\omega/\omega_{\rm crit}$ = {omega_over_omega_crit:.2f}')
        plt.title(r'Solution with optimal $j$ = ' +str(np.round(best_j, 2)) + ' for alpha = ' +  str(np.round(alpha_fixed,2)) + r' $\omega/\omega_{\rm crit}$ = ' + str(np.round( omega_over_omega_crit,2)))

        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        #plt.xlim(1, 2)
        
        #plt.xlim(.8,1)

        #plt.xlim(1.35, 1.6)
        #plt.ylim(1.35, 1.6)

        plt.show()
        
        
        
        d_omega_num = np.gradient(omega_sol, r)
        d_z_num = np.gradient(z_sol, r)
        f1_vals = np.array([f1(r[i], omega_sol[i], z_sol[i], best_j, alpha_fixed) for i in range(N)])
        f2_vals = np.array([f2(r[i], omega_sol[i], z_sol[i], best_j, alpha_fixed) for i in range(N)])
        
        
        
        d_omega_guess = np.gradient(omega_guess, r)
        d_z_guess = np.gradient(z_guess, r)
        f1_vals_guess = np.array([f1(r[i], omega_guess[i], z_guess[i], best_j, alpha_fixed) for i in range(N)])
        f2_vals_guess = np.array([f2(r[i], omega_guess[i], z_guess[i], best_j, alpha_fixed) for i in range(N)])
        
        
        
        
        # Plot derivative comparison
        #plt.subplot(1, 2, 1)
        plt.plot(r, d_omega_num, label=r'$\frac{d\omega}{dr}$ (numerical)', color='blue')
        plt.plot(r, f1_vals, '--', label=r'$f_1(r, \omega, z)$ (ODE RHS)', color='red')
        if True:
            plt.plot(r, d_omega_guess, label=r'$\frac{d\omega}{dr}$ (guess)', color='black')
            plt.plot(r, f1_vals_guess, '--', label=r'$f_1(r, \omega, z)$ (guess)', color='yellow')
        plt.xlabel('r')
        plt.title('Check: $\omega$ derivative')
        plt.legend()
        plt.ylim(-5, 10)
        
        plt.grid(True)
        
        
        
        fig, ax1=plt.subplots()        

        #plt.subplot(1, 2, 2)
        plt.plot(r, d_z_num, label=r'$\frac{dz}{dr}$ (numerical)', color='orange')
        plt.plot(r, f2_vals, '--', label=r'$f_2(r, \omega, z)$ (ODE RHS)', color='green')
        
        if True:
            plt.plot(r, d_z_guess, label=r'$\frac{dz}{dr}$ (guess)', color='black')
            plt.plot(r, f2_vals_guess, '--', label=r'$f_2(r, \omega, z)$ (guess)', color='yellow')
        plt.xlabel('r')
        plt.title('Check: $z$ derivative')
        plt.legend()
        plt.grid(True)
        
        
        plt.ylim(-2, 2)
        plt.tight_layout()
        plt.show()
        
        
        
    else:
        print("❌ Final solve with best j failed:", sol.message)
else:
    print("❌ Failed to optimize j:", result.message)












#%%
# === Solve the system ===
j = 0.5

# Initial guess: linear interpolation between BCs
omega0, z0 = inner_bc(j)
omega1, z1 = outer_bc(j)

omega_guess = np.linspace(omega0, omega1, N)
z_guess = np.linspace(z0, z1, N)
U0 = np.concatenate([omega_guess, z_guess])

# Solve
sol = root(relaxation_residual, U0, args=(j, alpha_fixed), method='lm', tol=1e-6)

# === Plot results ===
if sol.success:
    omega_sol = sol.x[:N]
    z_sol = sol.x[N:]

    plt.plot(r, omega_sol, label=r'$\omega(r)$', color='blue')
    plt.plot(r, z_sol, label=r'$z(r)$', color='orange')
    
    plt.plot(r, omega_guess, 'k--', label='Initial Guess', color='blue')
    plt.plot(r, z_guess, 'k--', label='Initial Guess', color='orange')

    
    
    
    plt.xlabel('r')
    plt.title('Solution via Relaxation Method')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
else:
    print("❌ Solver failed:", sol.message)