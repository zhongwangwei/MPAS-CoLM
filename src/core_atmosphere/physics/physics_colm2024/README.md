# MPAS-Embedded CoLM2024

This directory contains the CoLM2024 land-surface physics used by
MPAS-Atmosphere when `config_lsm_scheme = 'sf_colm2024'`.

Build CoLM2024 through MPAS:

- CMake: configure with `-DMPAS_COLM2024=ON`.
- Make: build atmosphere with `COLM2024=true`.

Run with `config_num_soil_layers = 10`. CoLM keeps its native 10-layer soil
column and its patch/PFT state; MPAS cells map one-to-one to CoLM elements.

MPAS initializes MPI and passes its domain communicator to CoLM. CoLM does not
start a separate process pool or require extra ranks. The MPAS LSM driver owns a
duplicated communicator context that isolates CoLM river-routing, remapping, and
restart messages while using the same MPAS processes. Distributed CoLM code still
issues the collectives and point-to-point exchanges required by those algorithms,
but only on this borrowed MPAS-owned context; it never owns an MPI runtime.

Standalone CoLM run cases, forcing namelist examples, and evaluation scripts
are intentionally not kept in this embedded source tree.
