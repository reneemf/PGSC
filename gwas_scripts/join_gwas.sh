#!/usr/bin/env bash

#SBATCH -J join_gwas
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --partition=tier2q
#SBATCH --mail-type ALL
#SBATCH --mail-user=<your_email>
#SBATCH -o /home/<username>/slurm_outputs/join_gwas/join_gwas_%A.out
#SBATCH -e /home/<username>/slurm_outputs/join_gwas/join_gwas_%A.err

now=$(date)
# Record job start time
start_time=$(date +%s)

echo "Date run: $now"

cd "${PGSC_HOME}"
source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"

join_gwas(){	
	pheno=$1
	context=$2
	# set up directories
	input_dir="${GWAS}${pheno}/"
	output_dir="${JOINED_GWAS}${pheno}/"
	mkdir -p "$output_dir"

	output_file="${output_dir}chr1-22_${TRAIN_POP}_${pheno}_${context}.assoc.linear"
	
	if [ ! -f "$output_file" ]; then
		first_file=$(ls ${input_dir}*${context}*.linear | head -1) 
		echo "$first_file"
		file_end=$(echo "$first_file" | cut -d'.' -f2-)
		chr1_gwas="${input_dir}chr1_${TRAIN_POP}_${pheno}_${context}.${file_end}"
		
		echo "joining ${chr1_gwas} and all other chromosomes"
		# write header to output file
		head -1 "$chr1_gwas" > "$output_file"
		# append chr data to output file
		tail -n +2 -q "${input_dir}"chr{1..22}_"${TRAIN_POP}"_"${pheno}"_"${context}"."${file_end}" >> "$output_file"
	else
		echo "${pheno} join gwas already made"
	fi
}

join_interaction_gwas(){
	pheno=$1
	context=$2
	context_code=$3
	# Create output directory if it doesn't exist
        input_dir="${GXEWAS}${pheno}/"
        scratch_output_dir="${SCRATCH_JOINED_GXEWAS}${pheno}/"
        output_dir="${JOINED_GXEWAS}${pheno}/"
	mkdir -p "$scratch_output_dir"
	mkdir -p "$output_dir"
	
        # Define file paths
	output_file="${scratch_output_dir}chr1-22_${TRAIN_POP}_${pheno}_${context}_GxC.assoc.linear"
	first_file=$(ls ${input_dir}*${context}_GxC*.linear | head -1)
        file_end=$(echo "$first_file" | cut -d'.' -f2-)
	chr1_gwas="${input_dir}chr1_${TRAIN_POP}_${pheno}_${context}_GxC.${file_end}"
	new_out_file="${output_dir}chr1-22_${TRAIN_POP}_${pheno}_${context}_GxC_crop.assoc.linear"

        if [ ! -f "$output_file" ]; then
		echo "joining ${chr1_gwas} and all other chromosomes"
		head -1 "$chr1_gwas" > "$output_file"
        tail -n +2 -q "${input_dir}"chr{1..22}_"${TRAIN_POP}"_"${pheno}"_"${context}"_GxC."${file_end}" >> "$output_file"
	else
		echo "original ${pheno} join gxewas already made"
	fi
	# extract header & keep only ADDxContext gwas output
	if [ ! -f "$new_out_file" ]; then
		head -1 "$output_file" > "$new_out_file"
		grep "ADDx${context_code}" "$output_file" >> "$new_out_file"
	else
		echo "cropped ${pheno} join gxewas already made"
	fi
}

context_code="" # code is the context's col name
for c in "${CONTEXTARRAY[@]}"
do
  echo "context: $c"
  # Process GWAS data for each phenotype
  for pheno in "${PHENOARRAY[@]}"; do
    if [[ "$c" = "sex" ]]; then
      context_code="31-0.0"
    elif [[ "$c" = "age" ]]; then
      context_code="34-0.0"
    else
      context_code="statins"
    fi

    # process the selected pheno's name
    parse_pheno=$(PROCESS_PHENO_NAME "$pheno")
    phenoLower=${parse_pheno##* }
		
    join_gwas "$phenoLower" "$c"
    join_interaction_gwas "$phenoLower" "$c" "$context_code"

  done
done

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds"
