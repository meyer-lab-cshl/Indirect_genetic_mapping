#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

#script to get the simulations
#libraries
library(optparse)
library(VariantAnnotation)
library(scran)
library(scuttle)
library(splatter)
library(stats)
library(tidyverse)

#reading in variables 
#params = readRDS("params_smartseq_use.rds")
#

#gtf file from Release 49 (GRCh38.p14) from https://www.gencodegenes.org/human/

option_list <- list(
  make_option("--vcf", type="character"),
  make_option("--kinmat", type="character"),
  make_option("--kinids", type="character"),
  make_option("--effect_size", type="double"),
  make_option("--samp_n", type="integer"),
  make_option("--cell_count", type="integer"),
  make_option("--rep", type="integer"),
  make_option("--gff", type="character"),
  make_option("--eqtl_params", type="character"),
  make_option("--rand", type="integer")
)

opt <- parse_args(OptionParser(option_list=option_list))

vcf_name <- opt$vcf
kin_mat_name <- opt$kinmat
kin_ids_name <- opt$kinids
effect_size <- opt$effect_size
samp_n <- opt$samp_n
cell_count <- opt$cell_count
rep <- opt$rep
gff <- read.table(opt$gff, sep = "\t")
eqtl_params <- readRDS(opt$eqtl_params)
set.seed(opt$rand)

path2plink <- "plink2"
path2ancids <- "../RESHAPE_geno_sim/ids/"
outputdir <- paste0("intermed_data/", "sim_geno", as.character(effect_size),"_sample",
                    as.character(samp_n), "_cell", as.character(cell_count), "_rep",
                    as.character(rep), "_rand", as.character(opt$rand)) 
#using parameter changes from the splatpop paper simulations github repo
eqtl_params <- setParams(eqtl_params,  
                         similarity.scale = 1.5,
                         batchCells = cell_count,
                         batch.size = samp_n,
                         nCells.sample = TRUE, 
                         nCells.shape = 1.387517805, 
                         nCells.rate = 0.006467842,
                         eqtl.dist = 1000000,
                         eqtl.maf.min = 0.1,
                         eqtl.maf.max = 0.5,
                         eqtl.ES.shape = 2.053687 + effect_size,
                         eqtl.ES.rate = 9.358627) 

#need to read in vcf file 
vcf <- readVcf(vcf_name, "hg38") 
#rownames(vcf) <- make.unique(rownames(vcf))
vcf <- vcf[!duplicated(rowRanges(vcf)), ]
vcf <- vcf[!duplicated(rownames(vcf)), ]

#this makes the key 
sim.means.A <- splatPopSimulateMeans(
  vcf = vcf,
  gff = gff,
  params = eqtl_params
)

#this uses the key to simulate count data
sim.sc.A <- splatPopSimulateSC(
  params = eqtl_params,
  key = sim.means.A$key,
  sim.means = sim.means.A$means,
  sparsify = FALSE
)

keyfilename = paste0(outputdir,"/key.rds")
saveRDS(sim.means.A, file = keyfilename)
###writing the files in the correct format for limix qtl####
#####phenotype file#####

#one file for each cell type 
#cell type A
# normalize within that group only
sim.sc.A <- logNormCounts(sim.sc.A)

pseudobulk <- aggregateAcrossCells(
  sim.sc.A,
  ids = colData(sim.sc.A)[, c("Sample")],
  use.assay.type = "logcounts",
  statistics = "mean"
)


#sim.sc <- logNormCounts(sim.sc)  
#norm_counts <- logcounts(sim.sc)
#pseudobulk <- aggregateAcrossCells(sim.sc, 
#                                   ids = colData(sim.sc)[,c("Sample", "Group")],
#                                   use.assay.type = "logcounts",  
#                                   statistics = "mean")  
#pheno_dataframe <- cbind(data.frame("feature_id" = rownames(pseudobulk)),
#                         data.frame(
#                           assay(pseudobulk[, pseudobulk$Group == "Group1"]))
#                         )
pheno_dataframe <- cbind(data.frame("feature_id" = rownames(pseudobulk)),
                         data.frame(assay(pseudobulk)))
pheno_dataframe$feature_id <- gsub("\\.", "-", pheno_dataframe$feature_id)

write.table(pheno_dataframe, file = paste0(outputdir, "/pheno_cellA.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

###################
####cell type B####
###################
#this makes the key 
sim.means.B <- splatPopSimulateMeans(
  vcf = vcf,
  gff = gff,
  params = eqtl_params
)

#this uses the key to simulate count data
sim.sc.B <- splatPopSimulateSC(
  params = eqtl_params,
  key = sim.means.B$key,
  sim.means = sim.means.B$means,
  method = "groups",
  sparsify = FALSE
)

# normalize within that group only
sim.sc.B <- logNormCounts(sim.sc.B)


#first selecting the 25% of genes to adjust in 
genes_to_mod <- sample(rownames(sim.sc.B), size = round(length(sim.sc.B) * 0.25))
#indirect genes
write.table(genes_to_mod, file = paste0(outputdir, "/indirect_gene.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


#overall hypothesis A:
#for cell type B, will select the same 25% of genes for both cell types
key <- sim.means.A$key
key$geneID <- rownames(key)

#group1_eqtls <- key[,
#  c("geneID", "eSNP.ID", "eQTL.EffectSize")
#]

group1_eqtls <- na.omit(key[, c("geneID", "eSNP.ID", "eQTL.EffectSize")])
geno_mat <- geno(vcf)$GT

gt_to_dosage <- function(gt) {
  gt_clean <- gsub("\\|", "/", gt)
  sapply(strsplit(gt_clean, "/"), function(alleles) sum(as.integer(alleles)))
}

dosage_mat <- apply(geno_mat, 2, gt_to_dosage)
snps_present <- group1_eqtls$eSNP.ID[group1_eqtls$eSNP.ID %in% rownames(dosage_mat)]
carrier_mat  <- dosage_mat[snps_present, , drop = FALSE] > 0

effect_sizes <- group1_eqtls$eQTL.EffectSize
names(effect_sizes) <- group1_eqtls$eSNP.ID

all_samples <- colnames(dosage_mat)

#gives you a matrix of the samples and the average eQTL simulated per sample 
#(i.e.) based on the genotypes of each sample
avg_effect_by_sample_overall <- setNames(
  sapply(all_samples, function(s) {
    carried <- snps_present[carrier_mat[snps_present, s]]
    if (length(carried) == 0) return(0)
    dosages <- dosage_mat[carried, s]
    weighted_effects <- abs(effect_sizes[carried]) * dosages
    sum(weighted_effects)
  }),
  all_samples
)

pseudobulk <- aggregateAcrossCells(
  sim.sc.B,
  ids = colData(sim.sc.B)[, c("Sample")],
  use.assay.type = "logcounts",
  statistics = "mean"
)

agg_g2 <- cbind(data.frame("feature_id" = rownames(pseudobulk)),
                         data.frame(
                           assay(pseudobulk))
)
print(colnames(agg_g2)
#colnames(agg_g2) <- c("feature_id", all_samples)

agg_g2_modified_overall <- agg_g2
agg_g2_modified_overall[genes_to_mod, all_samples] <-
  agg_g2[genes_to_mod, all_samples] +
  # broadcast the per-sample effect vector across all selected genes (row-wise)
  matrix(
    rep(avg_effect_by_sample_overall[all_samples], each = length(genes_to_mod)),
    nrow = length(genes_to_mod),
    ncol = length(all_samples),
    dimnames = list(genes_to_mod, all_samples)
  )

write.table(agg_g2_modified_overall, file = paste0(outputdir, "/pheno_cellB_modelA.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#cell type B (affected by the mean of eQTLs in cell type A) - pathway hypothesis
threshold <- quantile(abs(group1_eqtls$eQTL.EffectSize), 0.95)
top_eqtls <- group1_eqtls[abs(group1_eqtls$eQTL.EffectSize) >= threshold, ]

# Randomly select one
selected_eqtl <- top_eqtls[sample(nrow(top_eqtls), 1), ]


# Randomly select 4 additional SNPs from group1_eqtls excluding the selected eQTL
remaining_eqtls <- group1_eqtls[group1_eqtls$eSNP.ID != selected_eqtl$eSNP.ID, ]
additional_snps <- remaining_eqtls[sample(nrow(remaining_eqtls), 4), ]
write.table(c(additional_snps$eSNP.ID, selected_eqtl$eSNP.ID), 
            file = paste0(outputdir, "/modelB_selectedsnps.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ── Get samples that carry this specific eQTL's SNP 
selected_snp <- selected_eqtl$eSNP.ID

if (!selected_snp %in% rownames(dosage_mat)) {
  stop(sprintf("SNP %s not found in dosage matrix", selected_snp))
}

# Samples with dosage > 0 carry the variant
carrier_mask <- dosage_mat[selected_snp, ] > 0
carrier_samples_selected <- names(carrier_mask[carrier_mask])


# ── avg_effect_by_sample: only carrier samples, value is the selected eQTL's effect size
avg_effect_by_sample <- setNames(
  abs(selected_eqtl$eQTL.EffectSize) * dosage_mat[selected_snp, carrier_samples_selected],
  carrier_samples_selected
)

agg_g2_modified_pathway <- agg_g2
agg_g2_modified_pathway[genes_to_mod, carrier_samples_selected] <-
  agg_g2[genes_to_mod, carrier_samples_selected] +
  matrix(
    rep(avg_effect_by_sample[carrier_samples_selected], each = length(genes_to_mod)),
    nrow = length(genes_to_mod),
    ncol = length(carrier_samples_selected),
    dimnames = list(genes_to_mod, carrier_samples_selected)
  )
write.table(agg_g2_modified_pathway, file = paste0(outputdir, "/pheno_cellB_modelB.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#file for sample mapping
samplemap <- data.frame(geno_names = all_samples,
                        pheno_names = all_samples)

write.table(samplemap, file = paste0(outputdir,"/sample_map.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)

#for genetics files
vcfnewname = paste0(outputdir, "/vcf_to_plink.vcf")
writeVcf(vcf, filename = vcfnewname, 
         index = FALSE)

sys_com <- paste0(path2plink, " --vcf ", vcfnewname, 
                  " --vcf-half-call m --make-bed --out ", outputdir, "/geno")
system(sys_com)

##kinship file
kin_mat <- read.table(kin_mat_name)
kin_ids <- read.table(kin_ids_name, header = FALSE,
                      col.names = c("FID", "IID"))

double_kin_mat <- kin_mat * 2

labeled_kin_mat <- cbind(kin_ids$IID,
                         double_kin_mat)
bothlabel_kin_mat <- rbind(c("sample_id", kin_ids$IID), labeled_kin_mat)

write.table(bothlabel_kin_mat, file = paste0(outputdir,"/kinship.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = FALSE)

#covariate 
#first calculating the PCs
# PCs <- prcomp(assay(pseudobulk), rank = 20) 
# covar_unamed <- as.data.frame(PCs$rotation)
# covar_named <- cbind("sample_id" = rownames(covar_unamed),
#                      covar_unamed)
#then adding the ancestry, each ancestry has its own column 
#its 0 if the sample is not part of that ancestry and 1 if it is

#covariate 
#first calculating the PCs
#PCs <- prcomp(assay(pseudobulk), rank = 20) 
#covar_unamed <- as.data.frame(PCs$rotation)
covar_named <- cbind("sample_id" = all_samples,
                     batch = rep(1, samp_n))


write.table(covar_named, file = paste0(outputdir,"/covar.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

