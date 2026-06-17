# Edit these paths for your environment (or export PGSC_HOME in your shell).
UKB_SHARE <- "/path/to/UKB/dir/"          # shared UKB data root
SCRATCH   <- "/scratch/<username>/"           # your scratch space
PGSC_HOME <- Sys.getenv("PGSC_HOME", unset = paste0(UKB_SHARE,"analysis_dir/PGSC/"))

# Derived paths, no need to edit.
PHENOS <- paste0(UKB_SHARE,"extracted_phenotypes/")
COVARS <- paste0(PHENOS,"covariates_sa40PC/")
JOINED_GWAS <- paste0(PGSC_HOME,"pgsc_output/joined_gwas/")
JOINED_GXEWAS <- paste0(PGSC_HOME,"pgsc_output/joined_gxewas/")
PGS <- paste0(PGSC_HOME,"pgsc_output/pgs_out/")
PGS_GXE <- paste0(PGSC_HOME,"pgsc_output/pgs_gxe_out/")
R2 <- paste0(PGSC_HOME,"pgsc_output/r2_out/")

# Constants
TRAIN_POP <- "whitebrit"
TEST_POP <- "whitebrit"
VALID_POPS <- c("white_euro","asn","afr")
CONTEXT_LIST <- c("sex", "age", "statins")
THRESH <- c("0.0000000001", "0.00000001", "0.000001", "0.0001", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5")
PGS_COLS <- c("FID","IID","ALLELE_CT","ALLELE_DOSAGE_SUM","SCORE1_AVG","SCORE1_SUM")
PHENO_DIGITS <- as.character(commandArgs(TRUE)[[1]])
PHENO_NO_DIGITS <- ifelse(grepl("EA4", PHENO_DIGITS), "EA4", ifelse(grepl("FEV1674206", PHENO_DIGITS), "FEV1", 
                      ifelse(grepl("Arm_fat-free_mass_left_0", PHENO_DIGITS), "Arm_fat-free_mass_left_0", 
                          ifelse(grepl("FEV1_FVC_ratio674206", PHENO_DIGITS), "FEV1_FVC_ratio", 
                            ifelse(grepl("HbA1c674178", PHENO_DIGITS), "HbA1c", 
                              ifelse(grepl("Creatinine_2", PHENO_DIGITS), "Creatinine_2", 
                                ifelse(grepl("Height_0", PHENO_DIGITS), "Height_0", 
                                  ifelse(grepl("SHBG_0", PHENO_DIGITS), "SHBG_0", 
                                  ifelse(grepl("Bilirubin_0", PHENO_DIGITS), "Bilirubin_0", 
                                    ifelse(grepl("Lipoprotein_a_0", PHENO_DIGITS), "Lipoprotein_a_0", 
                                      ifelse(grepl("Apolipoprotein_a_0", PHENO_DIGITS), "Apolipoprotein_a_0", 
                                        ifelse(grepl("Apolipoprotein_b_0", PHENO_DIGITS), "Apolipoprotein_b_0", 
                                          ifelse(grepl("Arm_fat-free_mass_left_0", PHENO_DIGITS), "Arm_fat-free_mass_left_0", 
                                            ifelse(grepl("HDL_0", PHENO_DIGITS), "HDL_0", 
                                              ifelse(grepl("LDL_0", PHENO_DIGITS), "LDL_0", 
                                                ifelse(grepl("HbA1c_0", PHENO_DIGITS), "HbA1c_0", 
                                                  ifelse(grepl("Testosterone_0", PHENO_DIGITS), "Testosterone_0", 
                                                    ifelse(grepl("IGF-1674178", PHENO_DIGITS), "IGF-1", 
                                                      gsub("[0-9]+$", "", PHENO_DIGITS)))))))))))))))))))
PHENO_LOWER <- tolower(PHENO_NO_DIGITS)
SA_10_PC <- c("FID", "IID", "X31.0.0", "X34.0.0", "X22009.0.1", "X22009.0.2", "X22009.0.3", "X22009.0.4", "X22009.0.5", "X22009.0.6", "X22009.0.7", "X22009.0.8", "X22009.0.9", "X22009.0.10")
GWAS_COLS <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF?","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
CONTEXT_MAP <- list(sex     = list(col = "X31.0.0",     keep = 1:14),
                    age     = list(col = "X34.0.0",     keep = 1:14),
                    statins = list(col = "statins",     keep = c(1:14, 45)))

# Input files
PHENO_FILE <- paste0(PHENOS, PHENO_NO_DIGITS, "/", PHENO_DIGITS, ".pheno")

# Functions
R2fxn <- function(x,y,X){
  X <- as.matrix(X)
  y <- resid( lm( y ~ X ,na.action = na.exclude) )
  if( length(x) == length(y) ){
    x <- resid( lm( x ~ X ,na.action = na.exclude) )
  } else {
    for( j in 1:ncol(x) ){
      xj <- x[,j]
      x[,j] <- resid( lm( xj ~ X ,na.action = na.exclude) )
    }
  }
  summary( lm( y ~ 1 + x ) )
}


