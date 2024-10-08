#!/bin/bash

#SBATCH --job-name=hetero_sim   
#SBATCH --output=hetero_sim.out 
#SBATCH --time=06:00:00       
#SBATCH --ntasks=1            
#SBATCH --mem-per-cpu=3G      

module load r/4.3.2
# Run the R script with passed parameters
Rscript hetero_sim.R
