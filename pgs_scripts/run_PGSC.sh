#!/usr/bin/env bash

#SBATCH -J run_PGS
#SBATCH --mem=10GB
#SBATCH --time=0:10:00
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH -o /home/<username>/slurm_outputs/run_PGS/run_PGS_%A_%a.out
#SBATCH -e /home/<username>/slurm_outputs/run_PGS/run_PGS_%A_%a.err

now=$(date)
# Record job start time
start_time=$(date +%s)
echo "Date run: $now"

# make file containing the different P-value thresholds for inclusion of SNPs in the PGS
if [ ! -f std_range_list ]; then
	echo "0.0000000001 0 0.0000000001" > std_range_list
	echo "0.00000001 0 0.00000001" >> std_range_list
	echo "0.000001 0 0.000001" >> std_range_list
	echo "0.0001 0 0.0001" >> std_range_list
	echo "0.001 0 0.001" >> std_range_list
	echo "0.005 0 0.005" >> std_range_list
	echo "0.01 0 0.01" >> std_range_list
	echo "0.05 0 0.05" >> std_range_list
	echo "0.1 0 0.1" >> std_range_list
	echo "0.5 0 0.5" >> std_range_list
	echo "P-value thresholds range list created"
else
	echo "P-value thresholds range list already exists"
fi

sbatch test_PGS_PGSC.sh
#sbatch validate_PGS_PGSC.sh

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds" 

