# RG-D DIS Radiative Corrections

Reproducible inclusive DIS radiative-correction workflow for CLAS12 Run Group
D (RG-D). The repository contains the minimum EXTERNALS calculation engine,
RG-D target material definitions, kinematic run plans, execution/validation
scripts, and reviewed correction-factor tables.

The correction applied to a measured DIS-electron yield is

```text
C_DIS = sigma_Born / sigma_radiated
N_e_corrected = N_e_measured * C_DIS
```

## Origin and attribution

This work adapts the EXTERNALS package developed for inclusive electron
scattering analyses. The original repository is:

- https://github.com/utfsm-eg2-data-analysis/EXTERNALS

Suman Shrestha's historical reference fork is:

- https://github.com/suman0725/EXTERNALS

The RG-D repository is standalone so the production workflow can remain small
and auditable. Scientific provenance is preserved here and in `NOTICE.md`.

## Contents

```text
*.f, *.inc, Makefile   minimal EXTERNALS Fortran engine
INP/                   RG-D master inputs
TARG/                  LD2, C1, C2, Cu and Sn material definitions
RUNPLAN/               kinematic points supplied to EXTERNALS
scripts/               controlled execution and table extraction
results/               generated outputs (ignored until reviewed)
```

Legacy EG2 examples, Coulomb studies, scratch notebooks, compiled binaries,
object files, backup files, and old generated outputs are intentionally not
copied into this repository.

## Build on JLab ifarm

```bash
cd /work/clas12/suman/RGD_SIDIS_Analysis/projects/corrections/dis_rc/RGD-DIS-Radiative-Corrections
source set_env.sh
make
```

## Run the current RG-D reference grid

```bash
python3 scripts/run_rgd.py --run-plan RUNPLAN/rgd_kin.inp
```

The runner executes LD2, C1, C2, Cu and Sn independently. It retains the full
program log, extracts only the 13-column EXTERNALS physics rows, validates
positive Born/radiated cross sections, and writes a combined factor table.

Generated tables are placed under `results/`. They are not automatically
physics-approved or committed.

## Production requirement

The current `RUNPLAN/rgd_kin.inp` reproduces the earlier reference grid. The
final multiplicity-ratio analysis requires a separately reviewed run plan for
the exact hybrid xB-Q2 bins. Before publishing factors, verify:

1. beam energy and kinematic-point/bin-averaging prescription;
2. LD2 cell and window radiation lengths;
3. C1/C2/Cu/Sn foil thickness and upstream material ordering;
4. every Born and radiated cross section is finite and positive;
5. complete target and hybrid-bin coverage;
6. correction-factor systematic uncertainty.

## References

- L. W. Mo and Y. S. Tsai, *Rev. Mod. Phys.* **41**, 205 (1969).
- Y. S. Tsai, SLAC-PUB-848 (1971).

