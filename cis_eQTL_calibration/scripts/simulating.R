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
  make_option("--vcf_kin", type="character"),
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

path2plink <- "plink"
path2ancids <- "../RESHAPE_geno_sim/ids/"
outputdir <- paste0("intermed_data/", "sim_geno", as.character(effect_size),"_sample",
  as.character(samp_n), "_cell", as.character(cell_count), "_rep",
  as.character(rep)) 
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
                         eqtl.ES.shape = 2.053687,
                         eqtl.ES.rate = 9.358627) 

#need to read in vcf file 
vcf <- readVcf(vcf_name, "hg38") 
#rownames(vcf) <- make.unique(rownames(vcf))
vcf <- vcf[!duplicated(rowRanges(vcf)), ]
vcf <- vcf[!duplicated(rownames(vcf)), ]

#this makes the key 
sim.means <- splatPopSimulateMeans(
  vcf = vcf,
  gff = gff,
  params = eqtl_params
)

#this uses the key to simulate count data
sim.sc <- splatPopSimulateSC(
  params = eqtl_params,
  key = sim.means$key,
  sim.means = sim.means$means,
  sparsify = FALSE
)
print("key made")
#saving the key
saveRDS(sim.means$key, file = paste0(outputdir,"/key.rds"))

###writing the files in the correct format for limix qtl####
#####phenotype file#####
#normalizing counts
sim.sc <- logNormCounts(sim.sc)  
norm_counts <- logcounts(sim.sc)

pseudobulk <- aggregateAcrossCells(sim.sc, 
                                   ids = colData(sim.sc)[,c("Sample")],
                                   use.assay.type = "logcounts",  
                                   statistics = "mean")  
pheno_dataframe <- cbind(data.frame("feature_id" = rownames(pseudobulk)),
                         data.frame(assay(pseudobulk)))

write.table(pheno_dataframe, file = paste0(outputdir, "/sim_phenotype.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#covariate file
#determining PC

#then needing to get ancestry 

#file for sample mapping
sample_names <- colnames(pseudobulk)
samplemap <- data.frame(geno_names = sample_names,
                        pheno_names = sample_names)

write.table(samplemap, file = outputdir+"/bulksim_sample_map.txt",
            quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)

#for genetics files
sys_com <- paste0(path2plink, " --vcf ", vcf_name, " --make-bed --out ", outputdir, "/geno")
system(sys_com)

###feature annotation file 
feature_anot <- data.frame("feature_id" = sprintf("gene_%04d", 1:nrow(gff_gene)),
                           "chromosome" = gff$V3,
                           "start" = gff$V4,
                           "end" = gff$V5)
write.table(feature_anot, file = paste0(outputdir,"annot.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

##kinship file

kin_mat <- read.table(kin_mat_name)
kin_ids <- read.table(kid_ids_name)

double_kin_mat <- kin_mat * 2

colnames(double_kin_mat) <- kin_ids$V1
cbind("sample_id" = kin_ids$V1,
      double_kin_mat)

write.table(double_kin_mat, file = paste0(outputdir,"kinship.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

#covariate 
#first calculating the PCs
PCs <- prcomp(assay(pseudobulk), rank = 20) 
covar_unamed <- as.data.frame(PCs$rotation)
covar_named <- cbind("sample_id" = rownames(covar_unamed),
                     covar_unamed)
#then adding the ancestry, each ancestry has its own column 
#its 0 if the sample is not part of that ancestry and 1 if it is

#reading in the names of the samples
afr_samples <- read.table(paste0(path2ancids,"afr_ids.txt"))
amr_samples <- read.table(paste0(path2ancids,"amr_ids.txt"))
csa_samples <- read.table(paste0(path2ancids, "csa_ids.txt"))
eas_samples <- read.table(paste0(path2ancids,"eas_ids.txt"))
eur_samples <- read.table(paste0(path2ancids,"eur_ids.txt"))

afr_ancestry <- as.numeric(c(covar_named$sample_id %in% afr_samples))
amr_ancestry <- as.numeric(c(covar_named$sample_id %in% amr_samples))
csa_ancestry <- as.numeric(c(covar_named$sample_id %in% csa_samples))
eas_ancestry <- as.numeric(c(covar_named$sample_id %in% eas_samples))
eur_ancestry <- as.numeric(c(covar_named$sample_id %in% eur_samples))

covar_file <- cbind(covar_named, 
                    "afr_ancestry" = afr_ancestry,
                    "amr_ancestry" = amr_ancestry,
                    "csa_ancestry" = csa_ancestry,
                    "eas_ancestry" = eas_ancestry,
                    "eur_ancestry" = eur_ancestry)
write.table(covar_file, file = paste0(outputdir,"covar.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)
