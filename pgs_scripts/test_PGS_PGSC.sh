#!/usr/bin/env bash

#SBATCH -J PGS_PGSC_test
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=81GB
#SBATCH --partition=tier2q
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH -o /home/<username>/slurm_outputs/run_PGS/PGS_PGSC_%A_%a.out
#SBATCH -e /home/<username>/slurm_outputs/run_PGS/PGS_PGSC_%A_%a.err
#SBATCH --array=1-48 

start_time=$(date +%s)
echo "Date run: $(date)"

source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"
pop_path="${JOINED_GENOS}${TEST_POP}/"
split_path="${JOINED_GENOS}${TEST_POP}_split/"
pop_file="ukb_chr1-22_${TEST_POP}"
pgs_path="${PGS}${selected_pheno}/"
pgs_gxe_path="${PGS_GXE}${selected_pheno}/"

mkdir -p \
  "${CLUMPED_GENOS}" \
  "${pgs_path}" \
  "${pgs_gxe_path}" \
  "${JOINED_GWAS}snp_pval/"

# Calculate the index for the selected phenotype
pheno="${PHENOARRAY[$((SLURM_ARRAY_TASK_ID-1))]:-}"
if [[ -z $pheno ]]; then
    echo "no pheno present" 
    exit 1
fi

# process the selected pheno's name
parse_pheno=$(PROCESS_PHENO_NAME "$pheno")
selected_pheno=${parse_pheno##* }

echo "Building std ${selected_pheno} PGS"
for context in "${CONTEXTARRAY[@]}"; do
  echo "context: ${context}"  
  # clump snps
  output_file="${CLUMPED_GENOS}${pop_file}_${selected_pheno}_${context}"
  gwas_file="${JOINED_GWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}_${context}.assoc.linear"
  valid_snp_file="${CLUMPED_GENOS}${pop_file}_${selected_pheno}_${context}.valid.snp"
  if [[ ! -f "$output_file".clumps ]]; then
    ${LOCAL_PLINK2} \
      --bfile "${pop_path}${pop_file}" \
      --clump-p1 1 \
      --clump-r2 0.1 \
      --clump-kb 250 \
      --clump "$gwas_file" \
      --out "$output_file"
    # grab col 3 for rsids
    awk 'BEGIN{OFS=FS} NR!=1{print $3}' ${output_file}.clumps > "$valid_snp_file"
  fi

  # build pgs
  pgs_out_file="${pgs_path}${selected_pheno}_${context}_pgs"
  last_thresh="${THRESH[${#THRESH[@]}-1]}" #grab ex. thresh

  # for the --score flag : 3 7 12 grabs snp ID, a1 allele, & effect size estimate
  # for --q-score-range : cols 3 & 15 of "$gwas_file" grabs snp IDs & pval cols
  if [[ ! -f "$pgs_out_file"."$last_thresh".sscore ]]; then
    ${LOCAL_PLINK2} \
      --bfile "${pop_path}${pop_file}" \
      --keep "${split_path}${pop_file}_test" \
      --score "$gwas_file" 3 7 12 header cols=+scoresums ignore-dup-ids \
      --q-score-range std_range_list "$gwas_file" 3 15 header \
      --extract "$valid_snp_file" \
      --out "$pgs_out_file"
  fi

  # build pgxcs
  gxewas_file="${JOINED_GXEWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}_${context}_GxC_crop.assoc.linear"
  pgxcs_out_file="${pgs_gxe_path}${selected_pheno}_${context}_pgxcs"
  if [[ ! -f "$pgxcs_out_file"."$last_thresh".sscore ]]; then
    ${LOCAL_PLINK2} \
      --bfile "${pop_path}${pop_file}" \
      --keep "${split_path}${pop_file}_test" \
      --score "$gxewas_file" 3 7 12 header cols=+scoresums ignore-dup-ids \
      --q-score-range std_range_list "$gxewas_file" 3 15 header \
      --extract "$valid_snp_file" \
      --out "$pgxcs_out_file"
  fi
		
done
echo "test PGS done"

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds" 

