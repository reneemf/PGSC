#!/bin/bash

# Edit these paths for your environment (or export PGSC_HOME in your shell).
UKB_SHARE="/gpfs/data/ukb-share/"                        # shared UKB data root
SCRATCH="/scratch/<username>/"                            # your scratch space
PGSC_HOME="${PGSC_HOME:-${UKB_SHARE}dahl/<user>/PGSC/}"   # pipeline working directory
LOCAL_PLINK2="/path/to/plink2"                            # PLINK2 v2.00a6LM AVX2 Intel executable

# Derived paths, no need to edit.
PHENOS="${UKB_SHARE}extracted_phenotypes/"
COVARS="${PHENOS}covariates_sa40PC/"
GENOS="${SCRATCH}genotypes/"
JOINED_GENOS="${UKB_SHARE}genotypes/pop_genos/"
CLUMPED_GENOS="${SCRATCH}pgsc_out_scratch/clumped_genos/"
GWAS="${SCRATCH}pgsc_out_scratch/gwas/"
JOINED_GWAS="${PGSC_HOME}pgsc_output/joined_gwas/"
SCRATCH_JOINED_GWAS="${SCRATCH}pgsc_out_scratch/joined_gwas/"
GXEWAS="${SCRATCH}pgsc_out_scratch/gxewas/"
JOINED_GXEWAS="${PGSC_HOME}pgsc_output/joined_gxewas/"
SCRATCH_JOINED_GXEWAS="${SCRATCH}pgsc_out_scratch/joined_gxewas/"
PGS="${PGSC_HOME}pgsc_output/pgs_out/"
PGS_GXE="${PGSC_HOME}pgsc_output/pgs_gxe_out/"
R2="${PGSC_HOME}pgsc_output/r2_out/"

# Constants
TRAIN_POP="whitebrit"
TEST_POP="whitebrit"
VALID_POPS=("white_euro" "afr" "asn")
CONTEXTARRAY=("sex" "age" "statins")
THRESH=(0.0000000001 0.00000001 0.000001 0.0001 0.001 0.005 0.01 0.05 0.1 0.5)
PHENOARRAY=("Alkaline_phosphatase674206" "Alanine_aminotrans674206" "Apolipoprotein_a674206"
"Apolipoprotein_a_0" "Apolipoprotein_b674206" "Apolipoprotein_b_0"
"Arm_fat-free_mass_left674178" "Arm_fat-free_mass_left_0" 
"Arm_fat-free_mass_avg" "Arm_fat-free_mass_right674178" 
"Aspartate_aminotrans674206" "Basophill_count674178" "Bilirubin674206" "Bilirubin_0"
"Birth_weight674178" "BMI674178" "Calcium674178" "Cholesterol674178"
"Creatinine674178" "Creatinine_urine674206" "Cystatin_c674206" "DiastolicBP_auto674178"
"Eosinophill_count674178" "FEV1_FVC_ratio674206"
"Gamma_glutamyltransferase674178" "HbA1c674178" "HbA1c_0"
"HDL674178" "HDL_0" "Heel_bone_mineral_density_Tscore674206" "Height674178"
"Height_0" "Hip_circumference674178" "IGF-1674178" "LDL674178" "LDL_0"
"Leukocyte_count674178" "Lipoprotein_a674206" "Lipoprotein_a_0" "Lymphocyte_count674206"
"Mean_corpuscular_volume674178" "Monocyte_count674206" "Neutrophill_count674206" 
"Phosphate674206" "Platelet_count674178" "Platelet_volume674206" "Protein674206"
"Pulse_rate674178" "RBC674178" "Right_hand_grip_strength674178"
"SHBG674178" "SHBG_0" "Sodium_urine674206" "SystolicBP_auto674178" "Testosterone674178"
"Testosterone_0" "Triglycerides674178" "Urate674178" "Urea674178" "Vitamin_D674178"
"Waist_circumference674178" "Whole_body_fat_mass674178" "WHRadjBMI_Zhu"
"Rheumatoid_factor674206" "Glucose674178" "Potassium_urine674206" )  # 66 phenos
MANTARRAY=("Cholesterol674178" "HbA1c674178" "HbA1c_0" "Lipoprotein_a674206" "Lipoprotein_a_0" "Apolipoprotein_a674206" "Apolipoprotein_a_0" "Apolipoprotein_b674206" "Apolipoprotein_b_0" "Arm_fat-free_mass_left674178" "Arm_fat-free_mass_left_0" "Arm_fat-free_mass_avg674178" "HDL674178" "HDL_0" "Height674178" "Height_0" "LDL674178" "LDL_0" "SHBG674178" "SHBG_0" "Testosterone674178" "Testosterone_0" "Urate674178" "WHRadjBMI_Zhu")
SA_PC_10=("FID" "IID" "31-0.0" "34-0.0" "22009-0.1" "22009-0.2" "22009-0.3" "22009-0.4" "22009-0.5" "22009-0.6" "22009-0.7" "22009-0.8" "22009-0.9" "22009-0.10")

# Input files
STD_COVARS="${COVARS}/covariates_sa40PC674178.pheno"

# Functions
PROCESS_PHENO_NAME() {
        local pheno_name=$1
        local phenoNoDigits=${pheno_name%%[0-9]*}
        local phenoLower=$(echo "$phenoNoDigits" | tr '[:upper:]' '[:lower:]')

        if [[ $pheno_name == *"FEV1674206"* ]]; then
                phenoNoDigits="FEV1"
                phenoLower="fev1"
        elif [[ $pheno_name == *"FEV1_FVC_ratio674206"* ]]; then
                phenoNoDigits="FEV1_FVC_ratio"
                phenoLower="fev1_fvc_ratio"
        elif [[ $pheno_name == *"IGF-1674178"* ]]; then
                phenoNoDigits="IGF-1"
                phenoLower="igf-1"
	      elif [[ $pheno_name == *"HbA1c674178"* ]]; then
		            phenoNoDigits="HbA1c"
                phenoLower="hba1c"
        elif [[ $pheno_name == *"Height_0"* ]]; then
                phenoNoDigits="Height_0"
                phenoLower="height_0"
        elif [[ $pheno_name == *"EA4"* ]]; then
                phenoNoDigits="EA4"
                phenoLower="ea4"
        elif [[ $pheno_name == *"SHBG_0"* ]]; then
                phenoNoDigits="SHBG_0"
                phenoLower="shbg_0"
        elif [[ $pheno_name == *"HbA1c_0"* ]]; then
                phenoNoDigits="HbA1c_0"
                phenoLower="hba1c_0"
        elif [[ $pheno_name == *"LDL_0"* ]]; then
                phenoNoDigits="LDL_0"
                phenoLower="ldl_0"
        elif [[ $pheno_name == *"HDL_0"* ]]; then
                phenoNoDigits="HDL_0"
                phenoLower="hdl_0"
        elif [[ $pheno_name == *"Apolipoprotein_a_0"* ]]; then
                phenoNoDigits="Apolipoprotein_a_0"
                phenoLower="apolipoprotein_a_0"
        elif [[ $pheno_name == *"Apolipoprotein_b_0"* ]]; then
                phenoNoDigits="Apolipoprotein_b_0"
                phenoLower="apolipoprotein_b_0"
        elif [[ $pheno_name == *"Lipoprotein_a_0"* ]]; then
                phenoNoDigits="Lipoprotein_a_0"
                phenoLower="lipoprotein_a_0"
        elif [[ $pheno_name == *"Bilirubin_0"* ]]; then
                phenoNoDigits="Bilirubin_0"
                phenoLower="bilirubin_0"
        elif [[ $pheno_name == *"Arm_fat-free_mass_left_0"* ]]; then
                phenoNoDigits="Arm_fat-free_mass_left_0"
                phenoLower="arm_fat-free_mass_left_0"
        elif [[ $pheno_name == *"Testosterone_0"* ]]; then
                phenoNoDigits="Testosterone_0"
                phenoLower="testosterone_0"
	fi

        echo "$phenoNoDigits $phenoLower" 
}


