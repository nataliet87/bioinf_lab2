#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nt1598@nyu.edu
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=00-12:00:00
#SBATCH --output=serial_test_%j.log#!/bin/bash

module add kallisto/0.44.0
kallisto index -i /gpfs/scratch/nt1598/Homo_sapiens/Transcriptome/hg38transcripts.idx /gpfs/scratch/nt1598/Homo_sapiens/Transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz

kallisto quant -i /gpfs/scratch/nt1598/Homo_sapiens/Transcriptome/hg38transcripts.idx -o /gpfs/scratch/nt1598/data_bioinf/lab2/kallisto/BU1 /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU1_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU1_2.fastq.gz
kallisto quant -i /gpfs/scratch/nt1598/Homo_sapiens/Transcriptome/hg38transcripts.idx -o /gpfs/scratch/nt1598/data_bioinf/lab2/kallisto/BU2 /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU2_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU2_2.fastq.gz
kallisto quant -i /gpfs/scratch/nt1598/Homo_sapiens/Transcriptome/hg38transcripts.idx -o /gpfs/scratch/nt1598/data_bioinf/lab2/kallisto/BU3 /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU3_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU3_2.fastq.gz
kallisto quant -i /gpfs/scratch/nt1598/Homo_sapiens/Transcriptome/hg38transcripts.idx -o /gpfs/scratch/nt1598/data_bioinf/lab2/kallisto/DS1 /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H1_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H1_2.fastq.gz
kallisto quant -i /gpfs/scratch/nt1598/Homo_sapiens/Transcriptome/hg38transcripts.idx -o /gpfs/scratch/nt1598/data_bioinf/lab2/kallisto/DS2 /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H2_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H2_2.fastq.gz 
kallisto quant -i /gpfs/scratch/nt1598/Homo_sapiens/Transcriptome/hg38transcripts.idx -o /gpfs/scratch/nt1598/data_bioinf/lab2/kallisto/DS3 /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H3_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H3_2.fastq.gz
