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
keyfilename = paste0(outputdir,"/key.rds")
saveRDS(sim.means$key, file = keyfilename)

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
pheno_dataframe$feature_id <- gsub("\\.", "-", pheno_dataframe$feature_id)

write.table(pheno_dataframe, file = paste0(outputdir, "/pheno.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#covariate file
#determining PC

#then needing to get ancestry 

#file for sample mapping
sample_names <- colnames(pseudobulk)
samplemap <- data.frame(geno_names = sample_names,
                        pheno_names = sample_names)

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


#colnames(double_kin_mat) <- kin_ids$IID
labeled_kin_mat <- cbind(kin_ids$IID,
      double_kin_mat)
bothlabel_kin_mat <- rbind(c("sample_id", kin_ids$IID), labeled_kin_mat)

write.table(bothlabel_kin_mat, file = paste0(outputdir,"/kinship.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = FALSE)

#covariate 
#first calculating the PCs
#PCs <- prcomp(assay(pseudobulk), rank = 20) 
#covar_unamed <- as.data.frame(PCs$rotation)
covar_named <- cbind("sample_id" = sample_names,
                     batch = rep(1, samp_n))
#then adding the ancestry, each ancestry has its own column 
#its 0 if the sample is not part of that ancestry and 1 if it is

#reading in the names of the samples

write.table(covar_named, file = paste0(outputdir,"/covar.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)
