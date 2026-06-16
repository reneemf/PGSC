#!/usr/bin/env bash

#SBATCH -J sig_loci
#SBATCH --mem=10GB
#SBATCH --time=01:30:00
#SBATCH --partition=tier1q
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH -o /home/<username>/slurm_outputs/sig_loci/sig_loci_%A_%a.out
#SBATCH -e /home/<username>/slurm_outputs/sig_loci/sig_loci_%A_%a.err

now=$(date)
# Record job start time
start_time=$(date +%s)
echo "Date run: $now"

source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"
pop_file="ukb_chr1-22_${TEST_POP}"
sig_loci="${PGSC_HOME}pgsc_output/sig_loci/"
mkdir -p ${JOINED_GXEWAS}snp_pval/
mkdir -p ${sig_loci}

for context in ${CONTEXTARRAY[@]}; do
  sig_loci_file="${sig_loci}sig_loci_${context}.txt" # final output file
  for pheno in ${PHENOARRAY[@]}; do
    
    # process the selected pheno's name
    parse_pheno=$(PROCESS_PHENO_NAME "$pheno")
    selected_pheno=${parse_pheno##* }
    
    # filepaths
    gwas_file="${JOINED_GWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}_${context}.assoc.linear"
    gxewas_file="${JOINED_GXEWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}_${context}_GxC_crop.assoc.linear"
    valid_snp_file="${CLUMPED_GENOS}${pop_file}_${selected_pheno}_${context}.valid.snp"
    snps_pvals="${JOINED_GXEWAS}snp_pval/SNP.pvalue_${selected_pheno}_${context}"
    echo "$valid_snp_file"
    if [[ ! -f "$valid_snp_file" ]]; then # skip pheno
      echo "clumped valid snp file doesn't exist yet for pheno: $pheno, context: $context"
    else # grab sig loci
      # grab snp ID & pval cols from GxC gwas 
      awk '{print $3,$15}' ${gxewas_file} > ${snps_pvals}
      # grab sig snps from GxC gwas that are also in std clump
      awk 'NR==FNR{a[$1]; next} FNR==1 || ($1 in a) && ($2 < 5e-8)' ${valid_snp_file} ${snps_pvals} > ${sig_loci}best_snps/best_GxC_snps_${selected_pheno}_${context}.txt
      # count number of sig loci
      result=$(tail -n +2 ${sig_loci}best_snps/best_GxC_snps_${selected_pheno}_${context}.txt | wc -l)
      echo "${pheno} ${result}" >> ${sig_loci_file}
    fi
  done
done

# Record job end time
end_time=$(date +%s)

# Calculate total runtime in seconds
runtime=$((end_time - start_time))

# echo job information
echo "Job Runtime: $runtime seconds" 





