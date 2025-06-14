# Objective

Implement the ideas from [Paczynski
1991](https://ui.adsabs.harvard.edu/abs/1991ApJ...370..597P/abstract)
in MESA using `run_binary_extras.f90` physically motivated alternative
to rotationally-limited MT, possibly accounting for disk wind and L2
losses.

The relation \dot{M}(\dot{J}) from the Paczynski solution is computed
on the side and tabulated for application in MESA. For now interpolate
with `tanh` between Keplerian value at surface of accretor (if
accretor is non-critical) and 0 (if accretor approaches critical
rotation). When to interpolate is controlled by `x_ctrl` parameters in
`inlist_binary`. The non-conservativeness of mass transfer comes from
L2 losses, which also lead to losses of angular momentum.

Using [MESA r24.08.1](https://docs.mesastar.org/en/24.08.1/) or later
(`dev` version of MESA contain routines for L2 mass loss from Lu et
al. 2023)


# References or background material

-   [Paczynski 1991](https://ui.adsabs.harvard.edu/abs/1991ApJ...370..597P/abstract)
-   [Popham & Narayan 1991](https://ui.adsabs.harvard.edu/abs/1991ApJ...370..604P/abstract)
-   [Bisnovatyi-Kogan 1994](https://ui.adsabs.harvard.edu/abs/1994MNRAS.269..557B/abstract)
-   [Colpi et al. 1991](https://ui.adsabs.harvard.edu/abs/1991MNRAS.253...55C/abstract)
-   [Lu et al. 2022 (L2 mass loss)](https://academic.oup.com/mnras/article/519/1/1409/6886566)
-   [Ichikawa & Osaki, 1994](https://ui.adsabs.harvard.edu/abs/1994PASJ...46..621I/abstract)
-   [Papaloizou & Pringle 1977](https://academic.oup.com/mnras/article/181/3/441/988438)
-   [Paczynski 1977](https://ui.adsabs.harvard.edu/abs/1977ApJ...216..822P/abstract)
-   [Wave mediated angular momentum transport in astrophysical boundary layers](https://www.aanda.org/articles/aa/full_html/2015/07/aa26005-15/aa26005-15.html)
-   [On the terminal spins of accreting stars and planets: boundary layers](https://academic.oup.com/mnras/article/508/2/1842/6373455)
-   [The Effects of Cooling on Boundary Layer Accretion](https://arxiv.org/abs/2405.20367v1)
-   [Radiation hydrodynamics of the boundary layer in accretion disks. I - Numerical methods](https://ui.adsabs.harvard.edu/abs/1989A%26A...208...98K/abstract)
-   [Staristin 2022](https://ui.adsabs.harvard.edu/abs/2022RAA....22j5015S/abstract)

# Initial discussion:

-   the disk always fits within the Roche lobe (no truncation needed)
-   L2 mass loss is energetically possible (cf. [L2 mass loss
-   code](https://github.com/wenbinlu/L2massloss))

## Input & outputs

Input: $\dot{M}_1$, $\omega_{\rm accretor,surf}$, $M_1$, $M_2$, $a$

Output: $\dot{M}_2$, $\dot{J}_2$

## Routines to be used in `MESA`:

-   `mod_other_accreted_material_j.f90`
-   `mod_other_adjust_mdots.f90`
-   `mod_other_binary_jdot.f90`

see also `MESA_setup/src/binary_disk.inc` and `MESA_setup/src/run_binary_extras.f90`.

## Branches

- `default` contains a run without any custom implementation of L2 mass loss nor angular momentum accretion