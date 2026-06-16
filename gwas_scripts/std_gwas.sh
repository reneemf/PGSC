#!/usr/bin/env bash

#SBATCH -J std_gwas
#SBATCH --mem=20GB
#SBATCH --time=12:00:00
#SBATCH --partition=tier2q
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH --array=1-22%3
#SBATCH -o /home/<username>/slurm_outputs/std_gwas/std_gwas_%a_%A.out
#SBATCH -e /home/<username>/slurm_outputs/std_gwas/std_gwas_%a_%A.err

now=$(date)
# Record job start time
start_time=$(date +%s)

echo "Date run: $now"

source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"

cd ${PGSC_HOME}gwas_scripts/

pheno=$1
context=$2 
chr=${SLURM_ARRAY_TASK_ID}
pop_path="${GENOS}${TRAIN_POP}/"
pop_file="ukb_chr${chr}_${TRAIN_POP}_QC"
keep_file="${JOINED_GENOS}${TRAIN_POP}_split/ukb_chr1-22_${TRAIN_POP}_train_${context}"
echo "Processing pheno ${pheno} and chr ${chr}"

# process the selected pheno's name
parse_pheno=$(PROCESS_PHENO_NAME "$pheno")
phenoNoDigits=${parse_pheno%% *}
phenoLower=${parse_pheno##* }

# get pheno code
pheno_file="${PHENOS}${phenoNoDigits}/${pheno}.pheno"
pheno_code=$(head -1 "$pheno_file" | awk '{ print $3}')
gwas_path="${GWAS}${phenoLower}/"
mkdir -p ${gwas_path}

# pick cols given context type
if [[ "$context" = "sex" || "$context" = "age" ]]; then
  covar_file="${COVARS}covariates_sa40PC_${context}.pheno"
  ncovar='3-14'
else
  # context is statins
  covar_file="${COVARS}covariates_sa40PC_${context}.pheno"
  ncovar='3-14,45'
fi

# Define the output file path
out_file="${gwas_path}chr${chr}_${TRAIN_POP}_${phenoLower}_${context}"
if [[ ! -f "$out_file.$pheno_code.glm.linear" ]]; then
	${LOCAL_PLINK2} --bfile ${pop_path}${pop_file} \
	    --keep ${keep_file} \
	    --pheno "$pheno_file" \
	    --pheno-name "$pheno_code" \
	    --no-input-missing-phenotype \
	    --glm hide-covar \
	    --covar-variance-standardize 34-0.0 \
	    --covar "$covar_file" \
	    --covar-col-nums "$ncovar" \
	    --out "$out_file"
else
	echo "Output file already exists: $out_file.$pheno_code.glm.linear"
fi

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds"
