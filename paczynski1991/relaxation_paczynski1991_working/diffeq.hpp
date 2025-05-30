#ifndef diffeq_h
#define diffeq_h

struct Diffeq {
    Doub var_a_inv;
    Doub var_omega_spin;
    Doub var_n;
    VecDoub var_x_grid;
    Int mpt;
    Doub xbar, dx;
    Diffeq(const Doub &in_a_inv, const Doub &in_omega_spin, const Doub &in_n, const VecDoub in_x_grid);
    void smatrix(const Int k, const Int k1, const Int k2, const Int jsf, const Int is1, const Int isf, VecInt_I &indexv, MatDoub_O &s, MatDoub_I &y);
};

#endif
