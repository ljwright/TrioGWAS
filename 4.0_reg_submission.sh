#!/bin/bash -l
#$ -N triogwas_step_4
#$ -l h_rt=5:00:00
#$ -l mem=16G
#$ -pe smp 1
#$ -t 1-4352
#$ -wd "/home/rmjllwr/Scratch/Projects/MCS Trio GWAS/Liam Repo"

# Different batch controller at UCL
batch=$SGE_TASK_ID # Make a more expressive name
pheno=${1} # Tell people they need to pass the phenotype to this script.
# Could make it automated to batch to the right number of partitions - a submission script that takes $partitions as an argument
# Setting the working directory rather than just using it will defined file paths is a problem where the batch system saves error and output files to the working directory
#Working directory
# outm="/cluster/projects/trio_gwas" # Have to set this explicitly
outm="/home/rmjllwr/Scratch/Projects/MCS Trio GWAS/Liam Repo" # Does this need loading explicitly?

#Load R
#module add R/4.0.0-foss-2020a  
# module add R/4.2.1-foss-2022a
module -f unload compilers mpi gcc-libs # Have to load R already. Different on different systems
module load r/recommended


#Load Plink
# module add  plink/1.90b6.2  

#Call script - will pass through the phenotype name in each submission
chmod +x "${outm}/4.0_unified_regression.sh" # Need to change permissions
"${outm}/4.0_unified_regression.sh" $batch $pheno
echo "this is batch number $temp"
echo "this is pheno $pheno"


