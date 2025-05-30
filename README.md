# Objective

Implement the ideas from [Paczynski 1991](https://ui.adsabs.harvard.edu/abs/1991ApJ...370..597P/abstract) in MESA using
`run_binary_extras.f90` physically motivated alternative to
rotationally-limited MT, possibly accounting for disk wind and L2
losses.

**N.B.:** we should also use the [Lubow & Shu 1975](https://ui.adsabs.harvard.edu/abs/1975ApJ...198..383L/abstract) results to use this
 **only if** a disk does form (not for direct impact)

The relation \dot{M}(\dot{J}) from the Paczynski solution is computed
on the side and tabulated for application in MESA. The
non-conservativeness comes from L2 losses.

Using [MESA r24.08.1](https://docs.mesastar.org/en/24.08.1/) or later (`dev` version of MESA contain routines for
L2 mass loss from Lu et al. 2022)


# References or background material

-   [Paczynski 1991](https://ui.adsabs.harvard.edu/abs/1991ApJ...370..597P/abstract)
-   [Popham & Narayan 1991](https://ui.adsabs.harvard.edu/abs/1991ApJ...370..604P/abstract)
-   [Bisnovatyi-Kogan 1994](https://ui.adsabs.harvard.edu/abs/1994MNRAS.269..557B/abstract)
-   [Colpi et al. 1991](https://ui.adsabs.harvard.edu/abs/1991MNRAS.253...55C/abstract)
-   [Lu et al. 2022 (L2 mass loss)](https://academic.oup.com/mnras/article/519/1/1409/6886566)
-   [Ichikawa & Osaki, 1994](https://ui.adsabs.harvard.edu/abs/1994PASJ...46..621I/abstract)
-   [Papaloizou & Pringle 1977](https://academic.oup.com/mnras/article/181/3/441/988438)
-   [Paczynski 1977](https://ui.adsabs.harvard.edu/abs/1977ApJ...216..822P/abstract)


# Initial discussion:

-   the disk always fits within the Roche lobe (no truncation needed)
-   L2 mass loss is energetically possible (cf.
    [L2 mass loss code](https://github.com/wenbinlu/L2massloss)), but can be added later
-   the important thing to implement are dot(M2) and dot(J2) through the boundary layers.


## Input & outputs

Input: $\dot{M}_1$, $\omega_{\rm accretor,surf}$, $M_1$, $M_2$, $a$

Output: $\dot{M}_2$, $\dot{J}_2$

Routines to be used:

-   `mod_other_accreted_material_j.f90`
-   `mod_other_adjust_mdots.f90`
