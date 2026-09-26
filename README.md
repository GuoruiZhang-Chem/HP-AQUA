# HP-AQUA

This repository contains the 2-body, 3-body, and 4-body codes and fitted
coefficients used in HP-AQUA.

## Contents

- `2body`: H3O+ - H2O
- `3body`: H3O+ - (H2O)2
- `4body`: H3O+ - (H2O)3

The 2-body folder also contains the H3O+ PES/DMS source files and
coefficients needed by the current 2-body implementation.

The local absolute paths for coefficient files were changed to
relative paths.

The code was used with Intel Fortran (`ifort`) with `-O -r8`.

For the H3O+ PES/DMS library:

    cd 2body/src/h3o_pesdms_lib
    mkdir -p modh
    make lib

## Reference

G. Zhang, R. Ma, J. M. Bowman, and Q. Yu,
"Many-body machine-learned potential for the hydrated proton at
CCSD(T) accuracy with fast analytical gradient."

