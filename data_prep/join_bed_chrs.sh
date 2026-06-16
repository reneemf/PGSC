#!/usr/bin/env bash

#SBATCH --job-name=join_bed
#SBATCH --time=03:00:00
#SBATCH --mem=35gb
#SBATCH --partition=tier3q
#SBATCH --output=/home/<username>/slurm_outputs/join_bed/join_bed_%A.out
#SBATCH --error=/home/<username>/slurm_outputs/join_bed/join_bed_%A.err

source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"
cd ${GENOS}

### Function to merge files
merge_files() {
	local pop_type="$1"
  local mergelist="${pop_type}/mergelist.txt"
  local merged_prefix="${JOINED_GENOS}${pop_type}/ukb_chr1-22_${pop_type}"
	if [[ -e "${merged_prefix}.bed" ]]; then
		echo "Merged files already exist: ${merged_prefix}.*"
  else
	  rm -f "$mergelist"

    for i in {1..22}; do
      echo "${pop_type}/ukb_chr${i}_${pop_type}_QC" >> "$mergelist"
	  done

		${LOCAL_PLINK2} --pmerge-list "$mergelist" bfile \
      --make-bed \
      --out "$merged_prefix"
  fi
}

### merge training pop files
merge_files "$TRAIN_POP"

### merge validation pop files
for valid_pop in "${VALID_POPS[@]}"
do
  merge_files ${valid_pop}
done
