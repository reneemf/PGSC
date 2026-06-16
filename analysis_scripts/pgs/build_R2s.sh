#!/usr/bin/env bash

#SBATCH -J scale_process_pgs
#SBATCH --mem=40GB
#SBATCH --time=55:00:00
#SBATCH --partition=tier2q
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH -o /home/<username>/slurm_outputs/process_pgs/scale_process_pgs_%a_%A.out
#SBATCH -e /home/<username>/slurm_outputs/process_pgs/scale_process_pgs_%a_%A.err
#SBATCH --array=1-48


now=$(date)
# Record job start time
start_time=$(date +%s)

echo "Date run: $now"

source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"

# load required modules
module load gcc/12.1.0
module load R/4.3.1

# Calculate the index for the selected phenotype
pheno_num=$((SLURM_ARRAY_TASK_ID-1))
pheno="${PHENOARRAY[pheno_num]}"

if [[ ! $pheno ]]
then
        echo "no pheno present"
        exit
fi

Rscript compile_test_pgs.R "$pheno" 
#Rscript compile_valid_pgs.R "$pheno" 

# Record job end time
end_time=$(date +%s)

# Calculate total runtime in seconds
runtime=$((end_time - start_time))

# echo job information
echo "Job Runtime: $runtime seconds" 




