#####Libraries#####
library(optparse)
library(qvalue)
library(tidyverse)
library(data.table)

#####Variable loading#####
option_list <- list(
  make_option("--geneids", type="character"),
  make_option("--results", type="character"),
  make_option("--effect_size", type="double"),
  make_option("--samp_n", type="integer"),
  make_option("--cell_count", type="integer"),
  make_option("--rep", type="integer"),
  make_option("--rand", type="integer"),
  make_option("--model", type="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

outputfile <- paste0("results_celltype_B_model", as.character(opt$model), "/", as.character(opt$effect_size), 
                     "_sample", as.character(opt$samp_n), 
                     "_cell", as.character(opt$cell_count), "_rep",
                     as.character(opt$rep), "_rand", as.character(opt$rand), 
                     "_cellB", "/qtl_stats.txt")

keyname <- paste0("intermed_data/", "sim_geno", as.character(opt$effect_size),"_sample",
              as.character(opt$samp_n), "_cell", as.character(opt$cell_count), "_rep",
              as.character(opt$rep), "_rand", as.character(opt$rand),
              "/indirect_gene.txt") 


pathway_effect_name <- paste0("intermed_data/", "sim_geno", as.character(opt$effect_size),"_sample",
              as.character(opt$samp_n), "_cell", as.character(opt$cell_count), "_rep",
              as.character(opt$rep), "_rand", as.character(opt$rand),
              "/modelB_selectedsnps_effectsize.txt")

key <- read.table(keyname, header = FALSE)
pathway_effect <- read.table(pathway_effect_name, header = FALSE)

results_name <- opt$results
results <- data.table::fread(results_name)
fil_results <- results %>%
  group_by(feature_id) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(feature_id, snp_id, p_value, 
         indirect_molA_p_value, indirect_molA_effect_size, 
         indirect_molB_p_value, indirect_molB_effect_size)

fil_results$modelA_adjpval <- qvalue(fil_results$indirect_molA_p_value)$qvalue
fil_results$modelB_adjpval <- qvalue(fil_results$indirect_molB_p_value)$qvalue

sig_genes_modelA <- fil_results %>%
  filter(modelA_adjpval <= 0.05) %>%
  pull(feature_id) 


sig_genes_modelB <- fil_results %>%
  filter(modelB_adjpval <= 0.05) %>%
  pull(feature_id) 

#sum_results$modelA_adjpval <- qvalue(results$modelA_pval)$qvalue
#sig_genes_modelA <- sum_results[sum_results$modelA_adjpval <= 0.05, ]
  
#sum_results$modelB_adjpval <- qvalue(results$modelB_pval)$qvalue
#sig_genes_modelB <- sum_results[sum_results$modelB_adjpval <= 0.05, ]
ground_truth_ids <- key$V1



if (opt$model == "B"){ 
#  effect_size = pathway_effect$V1
  effect_size = opt$effect_size
} else {
  effect_size = opt$effect_size
}

#output files 
output_table_modelA <- data.frame("sim_data_model" = opt$model,
                           "mol_score_model" = "A",
                           "effect_size" = effect_size,
                           "sample_size" = opt$samp_n,
                           "cell_count" = opt$cell_count,
                           "rep" = opt$rep,
                           "FN" = length(setdiff(ground_truth_ids, sig_genes_modelA)),
                           "FP" = length(setdiff(sig_genes_modelA, ground_truth_ids)),
                           "TN" = nrow(fil_results) - length(union(sig_genes_modelA, ground_truth_ids)),
                           "TP" = length(intersect(sig_genes_modelA, ground_truth_ids))
                             )

output_table_modelB <- data.frame("sim_data_model" = opt$model,
                                  "mol_score_model" = "B",
                                  "effect_size" = effect_size,
                                  "sample_size" = opt$samp_n,
                                  "cell_count" = opt$cell_count,
                                  "rep" = opt$rep,
                                  "FN" = length(setdiff(ground_truth_ids, sig_genes_modelB)),
                                  "FP" = length(setdiff(sig_genes_modelB, ground_truth_ids)),
                                  "TN" = nrow(fil_results) - length(union(sig_genes_modelB, ground_truth_ids)),
                                  "TP" = length(intersect(sig_genes_modelB, ground_truth_ids))
)

output_table <- rbind(output_table_modelA, output_table_modelB)
write.table(output_table, file = outputfile, sep = ",", 
            row.names = FALSE, col.names = TRUE,
            quote = FALSE)

