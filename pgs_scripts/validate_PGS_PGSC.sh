#!/usr/bin/env bash

#SBATCH -J val_PGS_PGSC
#SBATCH --time=40:00:00
#SBATCH --mem-per-cpu=85GB
#SBATCH --partition=tier2q
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH -o /home/<username>/slurm_outputs/run_val_PGS/val_PGS_PGSC_%A_%a.out
#SBATCH -e /home/<username>/slurm_outputs/run_val_PGS/val_PGS_PGSC_%A_%a.err
#SBATCH --array=1-48

start_time=$(date +%s)
echo "Date run: $(date)"

source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"
best_thresh="${PGS}best_thresh/"
pgs_path="${PGS}valid/"
pgs_gxe_path="${PGS_GXE}valid/"

mkdir -p \
  "${best_thresh}" \
  "${pgs_path}" \
  "${pgs_gxe_path}"

# Calculate the index for the selected phenotype
pheno="${PHENOARRAY[$((SLURM_ARRAY_TASK_ID-1))]:-}"
if [[ -z $pheno ]]; then
  echo "no pheno present" 
  exit 1
fi
# process the selected pheno's name
parse_pheno=$(PROCESS_PHENO_NAME "$pheno")
selected_pheno=${parse_pheno##* }

for valid_pop in "${VALID_POPS[@]}"; do
  pop_path="${JOINED_GENOS}${valid_pop}/"
  pop_file="ukb_chr1-22_${valid_pop}"

  for context in "${CONTEXTARRAY[@]}"; do
  echo "Building std ${selected_pheno} PGSs, where context is ${context}, for ${valid_pop}"

    # grab best thresh from test pgs
    processed_pgs="${R2}processed_${selected_pheno}_${TEST_POP}_${context}_PGS.txt"
    pgs_valid_thresh="${best_thresh}${selected_pheno}_${context}_thresh.txt"
    best_pgs_thresh=$(head -2 "$processed_pgs" | tail -1 | awk '{print $2}' | tr -d '"')
    echo "$best_pgs_thresh 0 $best_pgs_thresh" > "$pgs_valid_thresh"
    echo "$best_pgs_thresh"

    # grab best generalized PGSC thresh
    pgsC_valid_thresh="${best_thresh}${selected_pheno}_${context}_pgxcs_thresh.txt"
    best_pgsC_thresh=$(head -4 "$processed_pgs" | tail -1 | awk '{print $2}' | tr -d '"')
    echo "$best_pgsC_thresh 0 $best_pgsC_thresh" > "$pgsC_valid_thresh"
    echo "$best_pgsC_thresh"  
    
    # build pgs
    pgs_out_file="${pgs_path}${selected_pheno}_${context}_${valid_pop}_pgs_valid"
    gwas_file="${JOINED_GWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}_${context}.assoc.linear"
    valid_snp_file="${CLUMPED_GENOS}ukb_chr1-22_${TEST_POP}_${selected_pheno}_${context}.valid.snp" #use test pgs clumped snps
  
    # pop filtering save one per phenoxcontext clumped snplist  
    if [[ ! -f ${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC1.bed ]]; then
      ${LOCAL_PLINK2} \
        --bfile "${pop_path}${pop_file}" \
        --extract "$valid_snp_file" \
        --make-bed \
        --out "${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC1"
    fi
    if [[ ! -f ${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC.bed ]]; then
      ${LOCAL_PLINK2} \
        --bfile "${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC1" \
        --mind 0.1 \
        --make-bed  \
        --out "${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC"
    fi
    # for the --score flag : 3 7 12 grabs snp ID, a1 allele, & effect size estimate
    # for --q-score-range : cols 3 & 15 of "$gwas_file" grabs snp IDs & pval cols
    if [[ ! -f "$pgs_out_file"."$best_pgs_thresh".sscore ]]; then
      ${LOCAL_PLINK2} \
        --bfile "${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC" \
        --keep "${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC.fam" \
        --score "$gwas_file" 3 7 12 header cols=+scoresums ignore-dup-ids \
        --q-score-range "$pgs_valid_thresh" "$gwas_file" 3 15 header \
        --extract "$valid_snp_file" \
        --out "$pgs_out_file"
    fi
    rm "${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC1"*

    # build gen pgxcs
    gxewas_file="${JOINED_GXEWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}_${context}_GxC_crop.assoc.linear"
    pgxcs_out_file="${pgs_gxe_path}${selected_pheno}_${context}_${valid_pop}_pgxcs_valid"
    if [[ ! -f "$pgxcs_out_file"."$best_pgsC_thresh".sscore ]]; then
      ${LOCAL_PLINK2} \
        --bfile "${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC" \
        --keep "${pop_path}QC/${pop_file}_${selected_pheno}_${context}_QC.fam" \
        --score "$gxewas_file" 3 7 12 header cols=+scoresums ignore-dup-ids \
        --q-score-range "$pgsC_valid_thresh" "$gxewas_file" 3 15 header \
        --extract "$valid_snp_file" \
        --out "$pgxcs_out_file"
    fi

  done
done

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds" 

