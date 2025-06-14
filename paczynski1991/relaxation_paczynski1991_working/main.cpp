#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>

#include "nr3.h"
#include "solvde.h"

#include "diffeq.hpp"

using namespace std;

Int main(Int argc, const char * argv[]) {
    Int N_grid = 1001;
    
    auto full_start = chrono::high_resolution_clock::now();
    float full_elapsed;
    
    VecDoub x_grid(N_grid);
    
    Doub h = (10.0 - 0.8) / ((Doub)N_grid-1.);
    for (Int ii=0; ii<N_grid; ii++) {
        x_grid[ii] = (Doub)ii * h + 0.8;
    }
    
    Int run_itmax = 1;
    Doub run_conv = 1e-15;
    Doub run_slowc = 0.;
    VecDoub run_scalv(3);
    run_scalv[0] = 10000.;
    run_scalv[1] = 10000.;
    run_scalv[2] = 1.;
    VecInt run_indexv(3);
    run_indexv[0] = 0;
    run_indexv[1] = 1;
    run_indexv[2] = 2;
    Int run_NB = 2;
    
    MatDoub y_vec(3,N_grid);    
    Int a_inv_points = 2;
    VecDoub a_inv_grid(a_inv_points);
    for (Int ii=0; ii<a_inv_points; ii++) {
        a_inv_grid[ii] = (2e-06 - 1e-06) * ((Doub)ii / ((Doub)a_inv_points - 1.)) + 1e-06;
    }
    
    Int omega_spin_points = 2;
    VecDoub omega_spin_grid(omega_spin_points);
    for (Int ii=0; ii<omega_spin_points; ii++) {
        omega_spin_grid[ii] = 0. * ((Doub)ii / ((Doub)omega_spin_points - 1.)) + 0.5; //0.54433105395;
    }
    
    Int n_points = 2;
    VecDoub n_grid(n_points);
    for (Int ii=0; ii<n_points; ii++) {
        n_grid[ii] = 1.5; //(1.5 - 1.5) * ((Doub)ii / ((Doub)n_points - 1.)) + 1.5;
    }
    
    ofstream outfile("/Users/nrui/Desktop/supercriticalmt/disk_RLOF_MESA/paczynski1991/relaxation_paczynski1991_working/output.txt");
    outfile.precision(10);
    outfile << "a_inv omega_spin n conv elapsed j_star";
    for (Int ii=0; ii<N_grid; ii++) {
        outfile << " omega" << to_string(ii);
    }
    for (Int ii=0; ii<N_grid; ii++) {
        outfile << " y" << to_string(ii);
    }
    outfile << '\n';
    outfile.close();
    
    Int convtest;
    float elapsed;
    float thresh = 0.001;
    
    for (Int a_inv_ii=0; a_inv_ii<a_inv_points; a_inv_ii++) {
        for (Int omega_spin_ii=0; omega_spin_ii<omega_spin_points; omega_spin_ii++) {
            for (Int n_ii=0; n_ii<n_points; n_ii++) {
                auto start = chrono::high_resolution_clock::now();
                convtest = 1;
                
                Diffeq diffeq(a_inv_grid[a_inv_ii], omega_spin_grid[omega_spin_ii], n_grid[n_ii], x_grid);
                
                for (Int ii=0; ii<N_grid; ii++) {
                }
                
                for (Int ii=0; ii<N_grid; ii++) {
                    y_vec[0][ii] = pow(x_grid[ii], -1.5);
                    y_vec[1][ii] = pow(a_inv_grid[a_inv_ii] / 1.5, 1./(2.*n_grid[n_ii]+3.)) * pow(x_grid[ii], (6.*n_grid[n_ii]+3.) / (4.*n_grid[n_ii]+6.));
                    y_vec[2][ii] = 0.1;
                }
                
                Solvde<Diffeq> run_solvde(run_itmax, run_conv, run_slowc, run_scalv, run_indexv, run_NB, y_vec, diffeq, convtest);
                
                auto stop = chrono::high_resolution_clock::now();
                elapsed = (float)((stop - start).count());                
                ofstream outfile("/Users/nrui/Desktop/supercriticalmt/disk_RLOF_MESA/paczynski1991/relaxation_paczynski1991_working/output.txt", std::ios_base::app);
                outfile.precision(10);
                
                outfile << to_string(a_inv_grid[a_inv_ii]) << " " << to_string(omega_spin_grid[omega_spin_ii]) << " " << to_string(n_grid[n_ii]) << " " << convtest << " " << elapsed << " " << to_string(y_vec[2][0]);
                for (Int jj=0; jj<N_grid; jj++) {
                    outfile << " " << to_string(y_vec[0][jj]);
                }
                for (Int jj=0; jj<N_grid; jj++) {
                    outfile << " " << to_string(y_vec[1][jj]);
                }
                outfile << '\n';
                outfile.close();
            }
        }
    }
    
    auto full_stop = chrono::high_resolution_clock::now();
    full_elapsed = (float)((full_stop - full_start).count());
    
    return 0;
}
