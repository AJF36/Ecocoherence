# Ecocoherence
This repository contains the scripts of the project ecological_coherence
Each of the folders contains the scripts for the function of the same name. For the scripts to work the structure of the repo must be respected

 
## Layout
- `analysis/scripts/` — pipeline scripts, numbered by dependency depth (`1_` … `8d_`) plus exploratory scripts. Run them in order; each reads its inputs from `analysis/results/<producer script>/`.
- `analysis/results/<script name>/`, `analysis/figures/<script name>/` — outputs of each script (git-ignored, regenerate with the scripts).
- `analysis/data/` — raw / source data (git-ignored).
- `analysis/reports/` — reports and manuscripts.
- `R/` — reusable functions sourced by the scripts.
- `scratch/` — work in progress, not tracked going forward.
- `docs/` — dependency diagram and per-script descriptions.

## Dependencies
- [fastspar](https://github.com/scwatts/fastspar/blob/main/README.md) for the creation of the correlation network
- [Functionink](https://github.com/apascualgarcia/functionInk) for the creation of modules. Functionink must be installed in the home/USER directory
- [anaconda](https://anaconda.org/) for managing the enviroments
- [carveme](https://carveme.readthedocs.io/en/latest/index.html) for building the metabolic models.
- [smetana](https://smetana.readthedocs.io/en/latest/) for calculating the cooperation/competition scores of the modules

