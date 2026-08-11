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
  make_option("--effect_size", type="double"),
  make_option("--samp_n", type="integer"),
  make_option("--cell_count", type="integer"),
  make_option("--rep", type="integer"),
  make_option("--eqtl_params", type="character"),
  make_option("--rand", type="integer")
)

opt <- parse_args(OptionParser(option_list=option_list))

kin_mat_name <- opt$kinmat
kin_ids_name <- opt$kinids
effect_size <- opt$effect_size
samp_n <- opt$samp_n
cell_count <- opt$cell_count
rep <- opt$rep
eqtl_params <- readRDS(opt$eqtl_params)
set.seed(opt$rand)

vcf <- mockVCF(
  n.snps = 100,
  n.samples = samp_n,
  chromosome = 1
)
gff <- mockGFF(
  n.genes = 50,
  chromosome = 1
)

path2plink <- "plink2"
outputdir <- paste0("intermed_data/", "sim_geno", as.character(effect_size),"_sample",
  as.character(samp_n), "_cell", as.character(cell_count), "_rep",
  as.character(rep), "_rand", as.character(opt$rand)) 


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
                         eqtl.ES.shape = 3.6 + effect_size,
                         eqtl.ES.rate = 12) 

sim.means <- splatPopSimulateMeans(
  vcf = vcf,
  gff = gff,
  params = eqtl_params
)
#single cell 
sim.sc <- splatPopSimulateSC(
  params = eqtl_params,
  key = sim.means$key,
  sim.means = sim.means$means,
  batchCells = cell_count,
  sparsify = FALSE
)

keyfilename = paste0(outputdir,"/key.rds")
saveRDS(sim.means$key, file = keyfilename)
print(keyfilename)
print(head(sim.means$key))
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

write.table(pheno_dataframe, file = paste0(outputdir, "/pheno.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#feature annot file
genomeanot <- gff
genomeanot$feature_id <- sprintf("gene_%02d", 1:50)
genomeanot <- genomeanot %>%
  mutate(chromosome = 1) %>%
  select(c(feature_id, chromosome, V4, V5))
colnames(genomeanot) <- c("feature_id", "chromosome",
                          "start", "end")
write.table(genomeanot, file = paste0(outputdir, "/genome_anot.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")



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

removing_old <- paste0("rm ", outputdir, "/vcf_to_plink.vcf")
system(removing_old)

##kinship file
fil_command = paste0("plink2 --bfile ", outputdir, "/geno --maf 0.05 --hwe 1e-6 --make-bed --out ",
                         outputdir, "/geno.fil")
prune_command = paste0("plink2 --bfile ", outputdir, "/geno.fil", " --indep-pairwise 250 50 0.2 --bad-ld --out ",
                       outputdir, "geno.fil.pruned")
kinship_mat = paste0("plink2 --bfile ", outputdir, "/geno.fil ", "--extract ",
                     outputdir, "geno.fil.pruned.prune.in  --make-king square --out ",
                     outputdir, "/king_ibd")
system(fil_command)
system(prune_command)
system(kinship_mat)

#putting it together in a kinship file
kin_mat <- read.table(paste0(outputdir, "/king_ibd.king"))
double_kin_mat <- kin_mat * 2
colnames(double_kin_mat) <- samplemap$geno_names
#double_kin_mat$sample_ids <- samplemap$sample_id

out_kin_mat <- cbind(sample_id = samplemap$geno_names, double_kin_mat)

write.table(out_kin_mat, file = paste0(outputdir,"/kinship.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

#covariate 
covar <- tibble(sample_id = colnames(vcf),
                batch = rep(1, samp_n))

write.table(covar, file = paste0(outputdir,"/covar.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)
