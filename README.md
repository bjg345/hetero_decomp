This repository contains the code used in the paper "A novel decomposition to explain heterogeneity in observational and randomized studies of causality." (Gilbert et al. 2024)

## Files in This Repository

- **`run.sh`**: A shell script used to run the main R script (`hetero_sim.R`). It uses the job scheduling system SLURM.

- **`submit_slurm.sh`**: A script that submits many jobs to SLURM for different parameter combinations. 

- **`hetero_sim.R`**: The main R script that runs the simulations. It calculates and breaks down heterogeneity between studies as described in the paper.

- **`functions.R`**: Contains helper functions used in the analysis.

- **`summarize.R`**: An R script that summarizes the simulation results and makes plots to show the findings.

