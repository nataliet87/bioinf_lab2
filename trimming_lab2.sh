#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nt1598@nyu.edu
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=00-12:00:00
#SBATCH --output=serial_test_%j.log#!/bin/bash

### Downloading runs from NCBI:
module add sratoolkit/2.9.1
fastq-dump --split-files SRR7049614 --gzip -O /gpfs/scratch/nt1598/data_bioinf

### loading needed packages for trimmming (trimgalore and cutadapt) and evaluating run quality:
module add miniconda3/4.5.1
conda activate /gpfs/scratch/nt1598/.conda/envs/seq_tools/
module add trimgalore/0.5.0
module add fastqc/0.11.7

# trims runs and generates fastqc reports:
trim_galore --illumina --pair --fastqc --output_dir /gpfs/scratch/nt1598/data_bioinf/lab2 ${1} ${2} &
trim_galore --illumina --pair --fastqc --output_dir /gpfs/scratch/nt1598/data_bioinf/lab2 ${3} ${4} &
trim_galore --illumina --pair --fastqc --output_dir /gpfs/scratch/nt1598/data_bioinf/lab2 ${5} ${6} &
trim_galore --illumina --pair --fastqc --output_dir /gpfs/scratch/nt1598/data_bioinf/lab2 ${7} ${8} &
trim_galore --illumina --pair --fastqc --output_dir /gpfs/scratch/nt1598/data_bioinf/lab2 ${9} ${10} &
trim_galore --illumina --pair --fastqc --output_dir /gpfs/scratch/nt1598/data_bioinf/lab2 ${11} ${12}
