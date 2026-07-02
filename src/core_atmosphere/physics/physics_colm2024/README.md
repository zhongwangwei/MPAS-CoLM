# MPAS-Embedded CoLM2024

This directory contains the CoLM2024 land-surface physics used by
MPAS-Atmosphere when `config_lsm_scheme = 'sf_colm2024'`.

Build CoLM2024 through MPAS:

- CMake: configure with `-DMPAS_COLM2024=ON`.
- Make: build atmosphere with `COLM2024=true`.

Standalone CoLM run cases, forcing namelist examples, and evaluation scripts
are intentionally not kept in this embedded source tree.
