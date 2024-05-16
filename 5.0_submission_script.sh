#!/bin/bash -l
#$ -m be
#$ -N triogwas_step_5
#$ -l h_rt=5:00:00
#$ -l mem=32G
#$ -pe smp 1
#$ -wd "/home/rmjllwr/Scratch/Projects/MCS Trio GWAS/Liam Repo"

phenotype=${1} # I think you could just write the script to loop over the available phenotypes
# I think this could be one script too as you aren't looping over anything.

#Working directory
outm="/home/rmjllwr/Scratch/Projects/MCS Trio GWAS/Liam Repo"

#Load R
module -f unload compilers mpi gcc-libs # Have to load R already. Different on different systems
module load r/recommended

#Call script 
chmod +x "${outm}"/5.0_combine_results.sh
"${outm}"/5.0_combine_results.sh $phenotype
