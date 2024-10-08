#!/bin/bash
#SBATCH --time=10:00:00         # Wall time limit (hh:mm:ss)

# Define parameter values
n_values=(100 1000 10000)
Q_values=(.1 .5)
B_values=(0 1)
C_values=(0 1)

# Loop through parameter values and submit jobs
for n in "${n_values[@]}"
do
  for Q_loop in "${Q_values[@]}"
  do
    Q_formatted=$(echo "$Q_loop" | sed 's/\.//') # Formatting Q value
    Q_formatted=$(printf "%02d" "$Q_formatted")
    
    for B_loop in "${B_values[@]}"
    do
      for C_loop in "${C_values[@]}"
      do
        for ((trial=1; trial<=1; trial++))
        do
          output_file="simulation_results/results_n${n}_Q${Q_formatted}_B${B_loop}_C${C_loop}_trial${trial}.rds"
          echo "Processing file: $output_file" 
          
          if [ -f "$output_file" ]; then
            echo "File $output_file already exists. Skipping..."
          else
            export E_N=$n
            export E_Q=$Q_loop
            export E_B=$B_loop
            export E_C=$C_loop
            export E_TRIAL=$trial

            sbatch run.sh

          fi
        done
      done
    done
  done
done
