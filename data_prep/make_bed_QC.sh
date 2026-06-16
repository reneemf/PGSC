#!/usr/bin/env bash 
#SBATCH --job-name=make_bed
#SBATCH --time=01:00:00
#SBATCH --mem=45gb
#SBATCH --partition=tier2q
#SBATCH --output=/home/<username>/slurm_outputs/make_bed/make_bed_%a_%A.out
#SBATCH --error=/home/<username>/slurm_outputs/make_bed/make_bed_%a_%A.err
#SBATCH --array=1-22

now=$(date)
# Record job start time
start_time=$(date +%s)

source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"
cd ${GENOS}

chr=${SLURM_ARRAY_TASK_ID} 

# ~/./plink2 is a local install, use plink2 if you don't have it locally but the flags may be different depending on version
${LOCAL_PLINK2} --pfile v3/ukb_imp_chr${chr}_v3 \
  --keep ${PHENOS}white_british/${TRAIN_POP}_unrelated.pheno \
  --make-bed \
  --out ${TRAIN_POP}/ukb_chr${chr}_${TRAIN_POP} 

# Delete SNPs & individuals w missingness, sig HWE pvals (genotyping error) & low MAF, remove duplicates & indels, name SNPs missing IDs (format: 'chr1:10583[b37]G,A')
# mind is a less stringent threshold to avoid being biased against Afr & Asn pops
${LOCAL_PLINK2} --bfile ${TRAIN_POP}/ukb_chr${chr}_${TRAIN_POP} \
  --maf 0.01 \
  --hwe 1e-6 \
  --geno 0.01 \
  --mind 0.1 \
  --rm-dup force-first \
  --set-missing-var-ids @:#[b37]\$r,\$a \
  --new-id-max-allele-len 50 missing \
  --make-bed  \
  --out ${TRAIN_POP}/ukb_chr${chr}_${TRAIN_POP}_QC

for valid_pop in "${VALID_POPS[@]}"
do
  mkdir -p ${valid_pop}
 
  # break up into valid pops, Remove duplicate SNPs, Name SNPs missing IDs (format: 'chr1:10583[b37]G,A')
  ${LOCAL_PLINK2} --pfile v3/ukb_imp_chr${chr}_v3 \
    --keep ${PHENOS}ethnic_background/${valid_pop}_unrelated.pheno \
    --rm-dup force-first \
    --set-missing-var-ids @:#[b37]\$r,\$a \
    --new-id-max-allele-len 50 missing \
    --make-bed \
    --out ${valid_pop}/ukb_chr${chr}_${valid_pop}_QC

done

# Record job end time
end_time=$(date +%s)

# Calculate total runtime in seconds
runtime=$((end_time - start_time))

# echo job information
echo "Job Runtime: $runtime seconds"
