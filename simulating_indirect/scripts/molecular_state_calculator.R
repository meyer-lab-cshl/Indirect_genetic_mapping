#script to determine the molecular score 

#####Libraries#####
library(optparse)
library(qvalue)
library(tidyverse)
library(data.table)
library(VariantAnnotation)
library(dplyr)

#####Variable loading#####
option_list <- list(
  make_option("--results", type="character"),
  make_option("--effect_size", type="double"),
  make_option("--samp_n", type="integer"),
  make_option("--cell_count", type="integer"),
  make_option("--rep", type="integer"),
  make_option("--rand", type="integer"),
  make_option("--threshold", type="double")
)

opt <- parse_args(OptionParser(option_list=option_list))

outputfile <- paste0("results_celltype_A/", as.character(opt$effect_size), "_sample",
                     as.character(opt$samp_n), "_cell", as.character(opt$cell_count), "_rep",
                     as.character(opt$rep), "_rand", as.character(opt$rand), "/qtl_stats.txt")
#naming variables
results_name <- opt$results
threshold <- opt$threshold
effect_size <- opt$effect_size
samp_n <- opt$samp_n
cell_count <- opt$cell_count
rep <- opt$rep

set.seed(opt$rand)

tempdir <- paste0("intermed_data/", "sim_geno", as.character(effect_size),"_sample",
                    as.character(samp_n), "_cell", as.character(cell_count), "_rep",
                    as.character(rep), "_rand", as.character(opt$rand)) 

vcfnewname = paste0(tempdir, "/vcf_to_plink.vcf")
vcf <- readVcf(vcfnewname, "hg38") 
#making the dosage matrix 
geno_mat <- geno(vcf)$GT

gt_to_dosage <- function(gt) {
  gt_clean <- gsub("\\|", "/", gt)
  sapply(strsplit(gt_clean, "/"), function(alleles) sum(as.integer(alleles)))
}

dosage_mat <- apply(geno_mat, 2, gt_to_dosage)
all_samples <- colnames(dosage_mat)

#getting eSNPs above the list 
results <- data.table::fread(results_name)
results <- na.omit(results)
results$qval <- qvalue(results$p_value)$qvalue
sig_results <- results[results$qval <= 0.05, ]

#model A - molecular score matrix
sig_snps_modelA <-  sig_results %>%
  dplyr::select(c(snp_id, beta)) %>%
  dplyr::filter(abs(beta) > threshold)

avg_beta_by_sample_modelA <- sapply(all_samples, function(s) {
  dosages <- dosage_mat[sig_snps_modelA$snp_id, s]
  sum(abs(sig_snps_modelA$beta) * dosages, na.rm = TRUE)
})

#when I test this out with a mock QTL, do the formatting for the output

#model B - 
#randomly subset for 5 eQTLs 
#(should I pick these from 5 randomly simulated eQTL maybe) - yes
# then of the five, see which ones are significant
#then average from there and print out 
selectedsnps_model_B_file <-  paste0(tempdir, "/modelB_selectedsnps.txt")

#will need to test out how to get the names of snps 
selectedsnps_model_B <- read.table(selectedsnps_model_B_file)
sig_snps_modelB <-  sig_results %>%
  dplyr::select(c(snp_id, beta)) %>%
  dplyr::filter(snp_id %in% selectedsnps_model_B$V1)

avg_beta_by_sample_modelB <- sapply(all_samples, function(s) {
  dosages <- dosage_mat[sig_snps_modelB$snp_id, s]
  sum(abs(sig_snps_modelB$beta) * dosages, na.rm = TRUE)
})

#again, when I test this out with a mock QTL, do the formatting for the output

covar_named <- cbind("sample_id" = names(avg_beta_by_sample_modelA),
                     "batch" = rep(1, samp_n),
                     "modelA" = unlist(unname(avg_beta_by_sample_modelA)),
                     "modelB" = unlist(unname(avg_beta_by_sample_modelB))
                       )

write.table(covar_named, file = paste0(tempdir,"/celltypeB_covar.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)
