#!/bin/env bash

#SBATCH --mem=40GB
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --job-name=nextflow_chip

#not using modules at rockefeller
#module load nextflow/23.04.4

source $HOME/.bashrc_rj_test.sh

conda activate nextflow_two
   ### not using this yet #SBATCH --mail-user=rjohnson@rockefeller.edu


# just checking if the pe_reads parameter can be seen form the workflow section.
# it can which means that this will work for any pair end reads for any project
# ues -entry parameter followed by the name of the workflow to run
#nextflow run chip_seq_nf_pipeline.nf  -resume --se_reads '../chip_fastqs/chip_*.fastq'

# i dont need to specify the bam files using the parameter since i hard coded them to be defualt
nextflow run nucleosome_nf_pipeline.nf -profile 'jan_peak_calling' -resume

