#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --mail-type=BEGIN,END,FAIL,STAGE_OUT
#SBATCH --mail-user=nt1598@nyu.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=00-12:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output=serial_test_%j.log#!/bin/bash

### loads version 2.2.9 of bowtie2 into my path (did this because I didn't realize I'd already downloaded an index, and encountered a bug
### when building index with newer and conda installed versions of bowtie2
export PATH=/gpfs/scratch/nt1598/packages/bowtie2-2.2.9:$PATH

### building index:
bowtie2-build -f --threads 2 /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa /gpfs/scratch/nt1598/Homo_sapiens/hgindex2/hg38idx

### aligning paired sequences against hg38 index; used preset sensitive parameters
bowtie2 --threads 2 --sensitive -x /gpfs/scratch/nt1598/Homo_sapiens/hgindex/hg38idx -1 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlBU1_1_val_1.fq.gz -2 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlBU1_2_val_2.fq.gz -S /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU1.sam &
bowtie2 --threads 2 --sensitive -x /gpfs/scratch/nt1598/Homo_sapiens/hgindex/hg38idx -1 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlBU2_1_val_1.fq.gz -2 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlBU2_2_val_2.fq.gz -S /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU2.sam &
bowtie2 --threads 2 --sensitive -x /gpfs/scratch/nt1598/Homo_sapiens/hgindex/hg38idx -1 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlBU3_1_val_1.fq.gz -2 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlBU3_2_val_2.fq.gz -S /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU3.sam &
bowtie2 --threads 2 --sensitive -x /gpfs/scratch/nt1598/Homo_sapiens/hgindex/hg38idx -1 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlDs12H1_1_val_1.fq.gz -2 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlDs12H1_2_val_2.fq.gz -S /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H1.sam &
bowtie2 --threads 2 --sensitive -x /gpfs/scratch/nt1598/Homo_sapiens/hgindex/hg38idx -1 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlDs12H2_1_val_1.fq.gz -2 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlDs12H2_2_val_2.fq.gz -S /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H2.sam &
bowtie2 --threads 2 --sensitive -x /gpfs/scratch/nt1598/Homo_sapiens/hgindex/hg38idx -1 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlDs12H3_1_val_1.fq.gz -2 /gpfs/scratch/nt1598/data_bioinf/lab2/CtrlDs12H3_2_val_2.fq.gz -S /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H3.sam
