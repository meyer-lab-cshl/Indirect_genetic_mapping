#!/bin/bash
#SBATCH -J filtering_gff
#SBATCH --mem=50G
#SBATCH --output=filtering_gff.o

Rscript filtering_gff.R
