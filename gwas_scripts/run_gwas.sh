now=$(date)
# Record job start time
start_time=$(date +%s)
echo "Date run: $now"
source "${PGSC_HOME:?Set PGSC_HOME to the repo directory}/source_file.sh"

# set up covariate files:
# update context values c_f = -(c_m*n_m)/n_f women=0, men=1 originally 
covar_sex="${COVARS}covariates_sa40PC_sex.pheno"
if [[ ! -f "$covar_sex" ]]; then
         # Count the number of 0s and 1s in column "31-0.0" (sex)
        n_f=$(awk '{if ($3 == 0) count++} END{print count}' "$STD_COVARS")
        n_m=$(awk '{if ($3 == 1) count++} END{print count}' "$STD_COVARS")
        # Calculate the replacement value
        replacement_value=$(awk -v N_m="$n_m" -v N_f="$n_f" 'BEGIN{printf "%.6f", -1*((1*N_m)/N_f)}')
        # Update the column and save the updated matrix
        awk -v x="$replacement_value" '{if ($3 == 0) $3 = x; print}' "$STD_COVARS" > "$covar_sex"
fi

covar_temp="${COVARS}covariates_sa40PC_temp.pheno"

covar_age="${COVARS}covariates_sa40PC_age.pheno"
if [[ ! -f "$covar_age" ]]; then
  # find mean yr of col 4 (UKB val: 1951.54): 1951.54
	avg=$(awk 'NR>1 {sum+=$4; count++} END {print sum/count}' "$STD_COVARS")
        # update values below and above mean in column "34-0.0" (age)
	head -n 1 "$STD_COVARS" > "$covar_temp"
	awk -v avg="$avg" 'NR>1 {if ($4 < avg) $4=0; else if ($4 > avg) $4=1; print}' "$STD_COVARS" >> "$covar_temp"
	# Count values below and above mean in column "34-0.0" (age)
	n_below=$(awk '{if ($4 == 0) count++} END{print count}' "$covar_temp")
        n_above=$(awk '{if ($4 == 1) count++} END{print count}' "$covar_temp")
	
        replacement_value=$(awk -v N_above="$n_above" -v N_below="$n_below" 'BEGIN{printf "%.6f", -1*((1*N_above)/N_below)}')
        
        awk -v x="$replacement_value" '{if ($4 == 0) $4 = x; print}' "$covar_temp" > "$covar_age"
	rm "$covar_temp"
fi

statins_pheno="${PHENOS}Statins/Statins674178.pheno"
covar_statins="${COVARS}covariates_sa40PC_statins.pheno"
if [[ ! -f "$covar_statins" ]]; then
        # add statin usage to covars file - it'll be col 45
        awk 'NR==FNR{a[$1]=$3; next} {if ($1 in a) print $0, a[$1]}' "$statins_pheno" "$STD_COVARS" > "$covar_temp"
        # Count the number of 0s and 1s in column "statins"
        n_no=$(awk '{if ($45 == 0) count++} END{print count}' "$covar_temp")
        n_yes=$(awk '{if ($45 == 1) count++} END{print count}' "$covar_temp")
       
        replacement_value=$(awk -v N_yes="$n_yes" -v N_no="$n_no" 'BEGIN{printf "%.6f", -1*((1*N_yes)/N_no)}')
        
        awk -v x="$replacement_value" '{if ($45 == 0) $45 = x; print}' "$covar_temp" > "$covar_statins"
        rm "$covar_temp"
fi

# actually push gwas & GxC gwas
for c in "${CONTEXTARRAY[@]}"
do
  # cut down genotypes
  keep_context_file="${JOINED_GENOS}${TRAIN_POP}_split/ukb_chr1-22_${TRAIN_POP}_train_${c}"
  if [[ ! -f "$keep_context_file" ]]; then
    covar_file="${COVARS}covariates_sa40PC_${c}.pheno"
    train_pop="${JOINED_GENOS}${TRAIN_POP}_split/ukb_chr1-22_${TRAIN_POP}_train"
    awk 'NR==FNR {ids[$1]; next} $1 in ids' "$covar_file" "$train_pop" > "$keep_context_file"
  fi

  for pheno in "${PHENOARRAY[@]}"; do
    echo "phenotype: $pheno"
    echo "context: $c"
    
    sbatch std_gwas.sh "$pheno" "$c"
    sleep 1
    sbatch gxe_gwas.sh "$pheno" "$c"
  done
done

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds"
