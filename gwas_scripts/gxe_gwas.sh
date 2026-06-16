#!/usr/bin/env bash

#SBATCH -J gxe_gwas
#SBATCH --mem=14GB
#SBATCH --time=08:00:00
#SBATCH --partition=tier2q
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH --array=1-22%3
#SBATCH -o /home/<username>/slurm_outputs/gxe_gwas/gxe_gwas_%a_%A.out
#SBATCH -e /home/<username>/slurm_outputs/gxe_gwas/gxe_gwas_%a_%A.err

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

# process the selected pheno's name
parse_pheno=$(PROCESS_PHENO_NAME "$pheno")
phenoNoDigits=${parse_pheno%% *}
phenoLower=${parse_pheno##* }

# get pheno code
pheno_file="${PHENOS}"${phenoNoDigits}"/"${pheno}".pheno"
pheno_code=$(head -1 "$pheno_file" | awk '{ print $3}')
gxewas_path="${GXEWAS}${phenoLower}/"
mkdir -p "${gxewas_path}"

# grab cols based on context type
echo "context: $context"
if [[ "$context" = "sex" ]]; then
  covar_file="${COVARS}covariates_sa40PC_${context}.pheno"
  ncovar='3-14'
  params='1-14'
elif [[ "$context" = "age" ]]; then
  covar_file="${COVARS}covariates_sa40PC_${context}.pheno"
  ncovar='3-14'
  params='1-13,15'  ##15 = addxage
else
  # context is statins
  covar_file="${COVARS}covariates_sa40PC_${context}.pheno"
  ncovar='3-14,45'
  params='1-14,27'  # 27 = addxstatins
fi

# Define the output file path
out_file="${gxewas_path}chr${chr}_${TRAIN_POP}_${phenoLower}_${context}_GxC"
if [[ ! -f "$out_file.$pheno_code.glm.linear" ]]; then
  ${LOCAL_PLINK2} --bfile ${pop_path}${pop_file} \
    --keep ${keep_file} \
    --pheno "$pheno_file" \
    --pheno-name "$pheno_code" \
    --no-input-missing-phenotype \
    --glm interaction \
    --parameters "$params" \
    --covar "$covar_file" \
    --covar-col-nums "$ncovar" \
    --covar-variance-standardize 34-0.0 \
    --out "$out_file"
else
  echo "Output file already exists: $out_file.$pheno_code.glm.linear"
fi	

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds"
