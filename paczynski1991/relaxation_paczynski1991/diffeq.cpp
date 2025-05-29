#include <cmath>
#include <vector>
#include "nr3.h"
#include "diffeq.hpp"

using namespace std;

Diffeq::Diffeq(const Doub &in_a_inv, const Doub &in_omega_spin, const Int &in_n, const VecDoub in_x_grid) {
    var_a_inv = in_a_inv;
    var_omega_spin = in_omega_spin;
    var_n = in_n;
    var_x_grid = in_x_grid;
    mpt = var_x_grid.size();
}

void Diffeq::smatrix(const Int k, const Int k1, const Int k2, const Int jsf, const Int is1, const Int isf, VecInt_I &indexv, MatDoub_O &s, MatDoub_I &y) {
    dx = var_x_grid[k] - var_x_grid[k-1];
    xbar = 0.5 * (var_x_grid[k] + var_x_grid[k-1]);
    
    if (k == k1) {
        s[1][jsf] = var_omega_spin - y[0][0];
        s[1][indexv[0]+3] = -1;
        s[1][jsf] = var_omega_spin - y[0][0];
        s[1][indexv[1]+3] = 0;
        s[1][jsf] = var_omega_spin - y[0][0];
        s[1][indexv[2]+3] = 0;
        s[2][jsf] = -y[1][0] + sqrt(pow(1 - 0.32000000000000006*pow(var_omega_spin, 2), -2.0) - 0.64000000000000012);
        s[2][indexv[0]+3] = 0;
        s[2][jsf] = -y[1][0] + sqrt(pow(1 - 0.32000000000000006*pow(var_omega_spin, 2), -2.0) - 0.64000000000000012);
        s[2][indexv[1]+3] = -1;
        s[2][jsf] = -y[1][0] + sqrt(pow(1 - 0.32000000000000006*pow(var_omega_spin, 2), -2.0) - 0.64000000000000012);
        s[2][indexv[2]+3] = 0;
        
    } else if (k > k2-1) {
        s[0][jsf] = 0.031622776601683791 - y[0][mpt-1];
        s[0][indexv[0]+3] = -1;
        s[0][jsf] = 0.031622776601683791 - y[0][mpt-1];
        s[0][indexv[1]+3] = 0;
        s[0][jsf] = 0.031622776601683791 - y[0][mpt-1];
        s[0][indexv[2]+3] = 0;
        
    } else {
        s[0][jsf] = var_a_inv*pow(xbar, 3.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0)*(-pow(xbar, 2)*((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k]) + (1.0/2.0)*y[2][k-1] + (1.0/2.0)*y[2][k]) - (-y[0][k-1] + y[0][k])/dx;
        s[0][indexv[0]] = -1.0/2.0*var_a_inv*pow(xbar, 2)*pow(xbar, 3.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0) + 1.0/dx;
        s[0][indexv[0]+3] = -1.0/2.0*var_a_inv*pow(xbar, 2)*pow(xbar, 3.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0) - 1/dx;
        s[0][jsf] = var_a_inv*pow(xbar, 3.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0)*(-pow(xbar, 2)*((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k]) + (1.0/2.0)*y[2][k-1] + (1.0/2.0)*y[2][k]) - (-y[0][k-1] + y[0][k])/dx;
        s[0][indexv[1]] = var_a_inv*pow(xbar, 3.0*var_n - 1.5)*(-1.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0)*(-pow(xbar, 2)*((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k]) + (1.0/2.0)*y[2][k-1] + (1.0/2.0)*y[2][k])/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]);
        s[0][indexv[1]+3] = var_a_inv*pow(xbar, 3.0*var_n - 1.5)*(-1.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0)*(-pow(xbar, 2)*((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k]) + (1.0/2.0)*y[2][k-1] + (1.0/2.0)*y[2][k])/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]);
        s[0][jsf] = var_a_inv*pow(xbar, 3.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0)*(-pow(xbar, 2)*((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k]) + (1.0/2.0)*y[2][k-1] + (1.0/2.0)*y[2][k]) - (-y[0][k-1] + y[0][k])/dx;
        s[0][indexv[2]] = (1.0/2.0)*var_a_inv*pow(xbar, 3.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0);
        s[0][indexv[2]+3] = (1.0/2.0)*var_a_inv*pow(xbar, 3.0*var_n - 1.5)*pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], -2.0*var_n - 3.0);
        s[1][jsf] = xbar*pow(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2), 1.5)*pow((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k], 2)/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]) - xbar/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]) - (-y[1][k-1] + y[1][k])/dx;
        s[1][indexv[0]] = xbar*pow(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2), 1.5)*((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k])/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]);
        s[1][indexv[0]+3] = xbar*pow(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2), 1.5)*((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k])/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]);
        s[1][jsf] = xbar*pow(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2), 1.5)*pow((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k], 2)/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]) - xbar/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]) - (-y[1][k-1] + y[1][k])/dx;
        s[1][indexv[1]] = xbar*sqrt(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2))*pow((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k], 2)*(0.75*y[1][k-1] + 0.75*y[1][k])/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]) - 1.0/2.0*xbar*pow(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2), 1.5)*pow((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k], 2)/pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2) + (1.0/2.0)*xbar/pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2) + 1.0/dx;
        s[1][indexv[1]+3] = xbar*sqrt(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2))*pow((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k], 2)*(0.75*y[1][k-1] + 0.75*y[1][k])/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]) - 1.0/2.0*xbar*pow(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2), 1.5)*pow((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k], 2)/pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2) + (1.0/2.0)*xbar/pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2) - 1/dx;
        s[1][jsf] = xbar*pow(pow(xbar, 2) + pow((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k], 2), 1.5)*pow((1.0/2.0)*y[0][k-1] + (1.0/2.0)*y[0][k], 2)/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]) - xbar/((1.0/2.0)*y[1][k-1] + (1.0/2.0)*y[1][k]) - (-y[1][k-1] + y[1][k])/dx;
        s[1][indexv[2]] = 0;
        s[1][indexv[2]+3] = 0;
        s[2][jsf] = (-y[2][k-1] + y[2][k])/dx;
        s[2][indexv[0]] = 0;
        s[2][indexv[0]+3] = 0;
        s[2][jsf] = (-y[2][k-1] + y[2][k])/dx;
        s[2][indexv[1]] = 0;
        s[2][indexv[1]+3] = 0;
        s[2][jsf] = (-y[2][k-1] + y[2][k])/dx;
        s[2][indexv[2]] = -1/dx;
        s[2][indexv[2]+3] = 1.0/dx;
    }
}
