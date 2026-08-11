library(GenomicRanges)
library(VariantAnnotation)

gff <- read.table("gencode.v49.annotation.chr2.gtf", sep = "\t")
vcf_name <- "intermed_data/sim_geno1.0_sample100_cell1000_rep3_rand466/snps_filtered_genome_chr2.vcf.gz"
gff_gene <- gff
gff_gene <- gff_gene[,1:8]
gff_gene <- gff_gene[gff_gene$V3=="gene",]

gene_gr <- GRanges(seqnames = gff_gene$V1, 
                   ranges = IRanges(start=gff_gene$V4, 
                                    end = gff_gene$V5))

#loading up variant info
vcf <- readVcf(vcf_name, "hg38") 
#vcf <- vcf[!isMultiAllelic(vcf)]
vcf <- vcf[!duplicated(rowRanges(vcf)), ]
vcf <- vcf[!duplicated(rownames(vcf)), ]
print("printing vcf numbers")
print(length(rownames(vcf)))
print(length(unique(rownames(vcf))))

variant_gr <- rowRanges(vcf)
variant_gr_100kb <- resize(
  variant_gr,
  width = width(variant_gr) + 10,
  fix = "center"
)

keep_genes <- overlapsAny(gene_gr, variant_gr_100kb)
gene_gr_filtered <- gene_gr[keep_genes]

#feature annotation file for limix
df_gff <- as.data.frame(gene_gr_filtered)

feat_anon <- data.frame("feature_id" = c(1:nrow(df_gff)),
                        "chromosome" = df_gff$seqname,
                        "start" = df_gff$start,
                        "end" = df_gff$end
                        )
write.table(feat_anon, file = "annot_chr2_onlygenesinrange.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

gff_file <- data.frame("V1" = df_gff$seqname,
                       "V2" = rep("source", nrow(df_gff)),
                       "V3" = rep("gene", nrow(df_gff)),
                       "V4" = df_gff$start,
                       "V5" = df_gff$end,
                       "V6" = rep(".", nrow(df_gff)),
                       "V7" = rep(".", nrow(df_gff)),
                       "V8" = rep(".", nrow(df_gff))
)

write.table(gff_file, file = "gff_chr2_onlygenesinrange.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

print("made gff files")

library(splatter)
eqtl_params <- readRDS("params_smartseq_use.rds")
eqtl_params <- setParams(eqtl_params,  
                         similarity.scale = 1.5,
                         batchCells = 1000,
                         batch.size = 100,
                         nCells.sample = TRUE, 
                         nCells.shape = 1.387517805, 
                         nCells.rate = 0.006467842,
                         eqtl.dist = 1000000,
                         eqtl.maf.min = 0.1,
                         eqtl.maf.max = 0.5,
                         eqtl.ES.shape = 2.053687,
                         eqtl.ES.rate = 9.358627) 
sim.means <- splatPopSimulateMeans(
  vcf = vcf,
  gff = gff_file,
  params = eqtl_params
)

sim.sc <- splatPopSimulateSC(
  params = eqtl_params,
  key = sim.means$key,
  sim.means = sim.means$means,
  sparsify = FALSE
)
print("sim done")
saveRDS(sim.means$key, file = paste0("test_key.rds"))
saveRDS(sim.sc, file = paste0("sc_sim.rds"))



