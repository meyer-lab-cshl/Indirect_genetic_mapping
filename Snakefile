import pandas as pd
import random
import sys

import os
from snakemake.io import glob_wildcards
random.seed(522)

effects = [0.75, 1, 1.25]
samp_numbers = [75, 100] #removed 125
cell_numbers = [50, 100, 500, 1000] #removing 1500, 250
replicates = [1, 2, 3]

ancestry = ["amr", "afr", "csa", "eas", "eur"]
gff_file = "gencode.v49.annotation.chr2.gtf"
eqtl_params ="params_smartseq_use.rds"

setup = pd.DataFrame(columns = ['effect_size', 'samp_n', 'cell_count', 'rep', 'rand_seed'])
for effect in effects:
    for sam in samp_numbers:
        for cell in cell_numbers:
            for replicate in replicates:
                rand = random.randint(1, 1000)
                row = {'effect_size':effect, 'samp_n':sam, 'cell_count':cell, 'rep':replicate, 'rand_seed':rand}
                row = pd.DataFrame([row])
                setup = pd.concat([setup, row])

rule all:
    input:
        "results/summary.tsv",
        expand("results/{setup.effect_size}_sample{setup.samp_n}_cell{setup.cell_count}_rep{setup.rep}_rand{setup.rand_seed}/qtl_stats.txt", 
             setup = setup.itertuples())
             
rule sample_ids:
    input:
        samples = "../RESHAPE_geno_sim/ids/{anc}_ids.txt"
    output: 
        randsamples = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/{anc}_sample_ids.txt"
    params:
        intermed = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}",
        num_per_anc = lambda wildcards: int(wildcards.samp_n) // 5
    resources:
        mem_mb = 1000
    shell:
        """
        mkdir -p {params.intermed}
        shuf {input.samples} \
        --random-source=<(python3 -c "import random, sys; random.seed({wildcards.rand_seed}); sys.stdout.buffer.write(random.randbytes(1000000))") \
        | head -n {params.num_per_anc} > {output.randsamples}
        """

rule filter_samples:
    input:
        sample_list = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/{anc}_sample_ids.txt",
        eQTLgenome = "../RESHAPE_geno_sim/sim/sim_{anc}_chr2.vcf.gz",
        kingenome = "../RESHAPE_geno_sim/sim/sim_{anc}_chr3.vcf.gz"
    output:
        eQTLchr = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/{anc}_genome_chr2.vcf.gz"),
        kinchr = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/{anc}_genome_chr3.vcf.gz")
    resources:
        mem_mb = 10000
    threads: 4
    shell:
        """
        bcftools view --threads {threads} -r chr2:1-125000000  -S {input.sample_list} {input.eQTLgenome} -Oz -o {output.eQTLchr}
        bcftools index -t {output.eQTLchr}

        bcftools view --threads {threads} -r chr3:1-10000000 -S {input.sample_list} {input.kingenome} -Oz -o {output.kinchr}
        bcftools index -t {output.kinchr}
        """

rule merge_genomes:
    input: 
        eqtlchr = lambda wc: expand(
            "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/{anc}_genome_chr2.vcf.gz",
            anc=ancestry,
            effect_size=wc.effect_size,
            samp_n=wc.samp_n,
            cell_count=wc.cell_count,
            rep=wc.rep,
            rand_seed=wc.rand_seed),
        kinchr = lambda wc: expand(
            "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/{anc}_genome_chr3.vcf.gz",
            anc=ancestry,
            effect_size=wc.effect_size,
            samp_n=wc.samp_n,
            cell_count=wc.cell_count,
            rep=wc.rep,
            rand_seed=wc.rand_seed
        )
    output:
        eqtlchr_merge = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/genome_chr2.vcf",
        kinchr_merge = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/genome_chr3.vcf")
    resources:
        mem_mb = 20000
    threads: 8
    shell:
        """
        bcftools merge --threads {threads} {input.eqtlchr} -O v -o {output.eqtlchr_merge}
        bcftools merge --threads {threads} {input.kinchr} -O v -o {output.kinchr_merge}
        """

rule filter_qc_genomes:
    input:
        eqtlchr_merge = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/genome_chr2.vcf",
        kinchr_merge = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/genome_chr3.vcf"
    output:
        eqtl_filter = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/filtered_genome_chr2.vcf.gz"),
        kin_filter = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/filtered_genome_chr3.vcf.gz")
    resources:
        mem_mb = 10000
    threads: 4
    shell:
        """
        vcftools --vcf {input.eqtlchr_merge} --maf 0.05 --hwe 0.001 --max-missing 0.95 --recode --recode-INFO-all --stdout | bgzip -c > {output.eqtl_filter}
        tabix -p vcf {output.eqtl_filter}
        vcftools --vcf {input.kinchr_merge} --maf 0.05 --hwe 0.001 --max-missing 0.95 --recode --recode-INFO-all --stdout | bgzip -c > {output.kin_filter}
        tabix -p vcf {output.kin_filter}
        """

rule sim_exp:
    input:
        sim_vcf = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/filtered_genome_chr2.vcf.gz",
        vcf_for_kin = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/filtered_genome_chr3.vcf.gz"
    output:
        pheno_file = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/pheno.txt"),
        covar = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/covar.txt"),
        kinship = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/kinship.txt"),
        sample_map = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/sample_map.txt"),
        key = temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/key.rds"),
        annotation= temp("intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/annot.txt")
    resources:
        mem_mb = 50000
    params:
        gff = gff_file,
        eqtl_params = eqtl_params
    log:
        "logs/sim_exp/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}.log"
    shell:
        """Rscript scripts/simulating.R \
            --vcf {input.sim_vcf} \
            --vcf_kin {input.vcf_for_kin} \
            --effect_size {wildcards.effect_size} \
            --samp_n {wildcards.samp_n} \
            --cell_count {wildcards.cell_count} \
            --rep {wildcards.rep} \
            --gff {params.gff} \
            --eqtl_param {params.eqtl_params} --rand {wildcards.rand_seed}"""

#snakemake rules from limix qtl 
## Settings
nGenes = 50 # recommend to set to 50; 5 used here for toy dataset
startPos = 0
endOffset = 1000000000
numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
windowSize = '1000000'
hwequilibrium = '0.000001'
FDR = '0.05'

rule generate_chunks:
    input:
        annotation="intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/annot.txt",
        pheno_file= "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/pheno.txt", 
    output:
        chunks="results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}" + "/chunks.txt"
    log:
        "logs/generate_chunks/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}.log"
    resources:
        mem_mb = 5000
    params:
        nGenes = nGenes,
        startPos = startPos,
        endOffset = endOffset
    script:
        "Limix_QTL/scripts/generate_chunks.R"

checkpoint setup_chunk_analysis:
    input:
        chunks= "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}" + "/chunks.txt"
    output:
        directory("results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}" + "/chunks")
    resources:
        mem_mb = 5000
    shell:
        """
        mkdir {output}
        for chunk in $(cat {input.chunks}); do
            echo $chunk > {output}/${{chunk//[:-]/_}}_info
        done
        """

rule run_qtl_mapping:
    input:
        chunkinfo = "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}" + "/chunks/{chunk}_info",
        af = gff_file,
        pf = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/pheno.txt",
        cf = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/covar.txt",
        kf = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/kinship.txt",
        smf = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/sample_map.txt",
        gen = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/geno.bed"
    output:
        "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}" + "/chunks/qtl_results_{chunk}.h5"
    log:
        "logs/run_qtl_mapping/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}_{chunk}.log"
    resources:
        mem_mb = 30000
    params:
        chunkstr=lambda wildcards: extendChunk(wildcards.chunk),
        od = "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}" + "/chunks",
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        hwe = hwequilibrium,
        w = windowSize,
    wildcard_constraints:
        chunk="\d{1,2}_\d*_\d*"
    shell:
        """
        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/run_QTL_analysis.py
        python Limix_QTL/run_QTL_analysis.py \
            --plink {input.gen} \
            -af {input.af} \
            -pf {input.pf} \
            -cf {input.cf} \
            -od {params.od} \
            -rf {input.kf} \
            --sample_mapping_file {input.smf} \
            -gr {params.chunkstr} \
            -np {params.np} \
            -maf {params.maf} \
            -hwe {params.hwe} \
            -w {params.w} \
            -c -gm gaussnorm -bs 500 -rs 0.95
        """

def collect_qtl_result_files(wildcards):
    checkpoint_output = checkpoints.setup_chunk_analysis.get(**wildcards).output[0]
    wc=glob_wildcards(os.path.join(checkpoint_output, "{chunk}_info"))
    return expand(checkpoint_output + "/qtl_results_{chunk}.h5",
            chunk=wc.chunk)

rule aggregate_qtl_results:
    input:
        collect_qtl_result_files,
    output:
        "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/qtl_results_all.txt"
    log:
        "logs/aggregate_qtl_results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}.log"
    resources:
        mem_mb = 10000
    params:
        IF = "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}" + "/chunks",
        OF = "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}" + "/"
    shell:
        """
        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py
        python Limix_QTL/post_processing/minimal_postprocess.py \
            -id {params.IF} \
            -od {params.OF} \
            -sfo
        """

rule determine_results:
    input:
        results = "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/qtl_results_all.txt",
        key = "intermed_data/sim_geno{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/key.rds"
    output:
        "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}/qtl_stats.txt"
    log:
        "logs/determine_results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}_rand{rand_seed}.log"
    resources:
        mem_mb = 3000
    shell:
        """
        Rscript scripts/summarizing_results.R \
            --key {input.key} \
            --results {input.results} \
            --effect_size {wildcards.effect_size} \
            --samp_n {wildcards.samp_n} \
            --cell_count {wildcards.cell_count} \
            --rep {wildcards.rep}
        """

rule agg_results:
    input:
        expand("results/{setup.effect_size}_sample{setup.samp_n}_cell{setup.cell_count}_rep{setup.rep}_rand{setup.rand_seed}/qtl_stats.txt",
        setup = setup.itertuples())
    output:
        "results/summary.tsv"
    resources:
        mem_mb = 3000
    log:
        "logs/agg_results/agg_results.log"
    shell:
        "scripts/agg_results.py" 


         
#old stuff - not needed
#rule top_feature:
#    input:
#        finalFile = "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}/qtl_results_all.txt"
#    output:
#        "results/{effect_size}_sample{samp_n}_cell{cell_count}_rep{rep}/top_qtl_results_all.txt"
#    params:
#        IF = "results/{wildcards.effect_size}_sample{wildcards.samp_n}_cell{wildcards.cell_count}_rep{wildcards.rep}" + "/chunks",
#        OF = "results/{wildcards.effect_size}_sample{wildcards.samp_n}_cell{wildcards.cell_count}_rep{wildcards.rep}" + "/"
#    shell:
#        """
#        #singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py
#        python Limix_QTL/post_processing/minimal_postprocess.py \
#            -id {params.IF} \
#            -od {params.OF} \
#            -tfb \
#            -sfo
#        """

 
