# count distribution of Cs across contexts & output table
source(file.path(Sys.getenv("PGSC_HOME"), "source_file.R"))

# loop through contexts
for(context in CONTEXT_LIST){
  # set up files
  genos <- paste0(UKB_SHARE,"genotypes/pop_genos/")
  covars <- read.table(paste0(COVARS,"covariates_sa40PC_",context,".pheno"),header = TRUE)
  file_map <- list(
    "White British train" = paste0(genos,"whitebrit_split/ukb_chr1-22_whitebrit_train_",context),
    "White British test" = paste0(genos,"whitebrit_split/ukb_chr1-22_whitebrit_test"),
    "European" = paste0(genos,"white_euro/ukb_chr1-22_white_euro.fam"),
    "African" = paste0(genos,"afr/ukb_chr1-22_afr.fam"),
    "Asian" = paste0(genos,"asn/ukb_chr1-22_asn.fam")
  )
  
  # grab context col from covar file & pick C types
  col_num <- 0
  C0_name <- ""
  C1_name <- ""
  if(context == "sex"){
    col_num <- 3
    C0_name <- "num female"
    C1_name <- "num male"
  } else if(context == "age"){
    col_num <- 4
    C0_name <- "num below mean age"
    C1_name <- "num above mean age"
  } else {
    col_num <- 45
    C0_name <- "num no statin use"
    C1_name <- "num statin user"
  }

  # Initialize table # to do: set up col names based on context type (maybe edit at end?)
  table_data <- data.frame(matrix(ncol = 3, nrow = 5))
  colnames(table_data) <- c("num individuals", C0_name,C1_name)
  rownames(table_data) <- names(file_map)
  
  # Fill in context counts based on matching IDs
  for (pop in names(file_map)) {
    # Read WB/EUR/AFR/ASN files
    ids <- read.table(file_map[[pop]], header = FALSE)[, 1]

    # Match IDs with covariates
    matched_ids <- covars[covars[, 1] %in% ids, ]
    print(head(covars[, col_num]))
    # Count
    num_indi <- dim(matched_ids)[1] 
    num_C1 <- sum(matched_ids[, col_num] == 1, na.rm = TRUE)
    num_C0 <- num_indi - num_C1

    # Update table
    table_data[pop, "num individuals"] <- num_indi
    table_data[pop, C0_name] <- num_C0
    table_data[pop, C1_name] <- num_C1
  }
  out_dir <- paste0(PGSC_HOME,"pgsc_output/pop_counts/")
  dir.create(out_dir, showWarnings = F)
  table_output <- paste0(out_dir,"pop_counts_",context,".txt")

  # Save table to file
  write.table(table_data, file = table_output)#, sep = "\t",quote = FALSE)
}



