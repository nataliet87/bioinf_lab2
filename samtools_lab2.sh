#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --mail-type=BEGIN,END,FAIL 
#SBATCH --mail-user=nt1598@nyu.edu
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=00-12:00:00
#SBATCH --output=serial_test_%j.log#!/bin/bash

module add samtools/1.9
module add bedtools/2.26.0

# parse sam files into zipped bam files
samtools view -b -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU1.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU1.sam & 
samtools view -b -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU2.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU2.sam &
samtools view -b -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU3.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU3.sam

samtools view -b -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H1.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H1.sam &
samtools view -b -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H2.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H2.sam &
samtools view -b -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H3.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H3.sam

# sort files
samtools sort -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU1.bam & 
samtools sort -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU2.bam &
samtools sort -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_BU3.bam
samtools sort -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H1.bam &
samtools sort -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H2.bam &
samtools sort -o /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/hg_align_Ds12H3.bam    

# index files
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU1out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU2out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU3out.bam
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H1out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H2out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H3out.bam

# extract reads from fragments mapping to foward strand

# first foward fragment:
samtools view -b -f99 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU1out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_BU1out.bam &
samtools view -b -f99 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU2out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_BU2out.bam &
samtools view -b -f99 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU3out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_BU3out.bam 
samtools view -b -f99 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H1out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_Ds12H1out.bam &
samtools view -b -f99 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H2out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_Ds12H2out.bam &
samtools view -b -f99 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H3out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_Ds12H3out.bam 

# second forward fragment:
samtools view -b -f147 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU1out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_BU1out.bam &
samtools view -b -f147 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU2out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_BU2out.bam &
samtools view -b -f147 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU3out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_BU3out.bam 
samtools view -b -f147 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H1out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_Ds12H1out.bam &
samtools view -b -f147 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H2out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_Ds12H2out.bam &
samtools view -b -f147 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H3out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_Ds12H3out.bam 
# merge first and second forward fragments:
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_BU1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_BU1out.bam &
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_BU2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_BU2out.bam &
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_BU3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_BU3out.bam
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_Ds12H1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_Ds12H1out.bam &
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_Ds12H2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_Ds12H2out.bam &
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd1_Ds12H3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/fwd2_Ds12H3out.bam

# extract reads from fragments mapping to reverse strand
# first reverse fragment:
samtools view -b -f83 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU1out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_BU1out.bam &
samtools view -b -f83 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU2out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_BU2out.bam &
samtools view -b -f83 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU3out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_BU3out.bam 
samtools view -b -f83 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H1out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_Ds12H1out.bam &
samtools view -b -f83 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H2out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_Ds12H2out.bam &
samtools view -b -f83 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H3out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_Ds12H3out.bam 

# second reverse fragment:
samtools view -b -f163 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU1out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_BU1out.bam &
samtools view -b -f163 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU2out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_BU2out.bam &
samtools view -b -f163 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU3out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_BU3out.bam 
samtools view -b -f163 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H1out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_Ds12H1out.bam &
samtools view -b -f163 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H2out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_Ds12H2out.bam &
samtools view -b -f163 /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H3out.bam > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_Ds12H3out.bam 
# merge first and second reverse fragments:
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_BU1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_BU1out.bam &
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_BU2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_BU2out.bam &
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_BU3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_BU3out.bam
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_Ds12H1out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_Ds12H1out.bam &
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_Ds12H2out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_Ds12H2out.bam &
samtools merge -f /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev1_Ds12H3out.bam /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/rev2_Ds12H3out.bam

# index extracted fragment files
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU1out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU2out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU3out.bam
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU1out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU2out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU3out.bam
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/foward_Ds12H1out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/foward_Ds12H2out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/foward_Ds12H3out.bam
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H1out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H2out.bam &
samtools index /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H3out.bam

# use bedtools to parse files into 4 column bedgraph files (Gviz can read these)
# BU1 files:
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU1out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU1.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU1out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU1.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU1out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/BU1.bedgraph 
# BU2 files:
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU2out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU2.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU2out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU2.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU2out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/BU2.bedgraph 
# BU3 files:
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU3out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_BU3.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU3out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_BU3.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_BU3out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/BU3.bedgraph 

# Ds12H1 files:
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H1out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H1.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H1out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H1.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H1out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/Ds12H1.bedgraph 
# Ds12H2 files:
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H2out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H2.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H2out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H2.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H2out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/Ds12H2.bedgraph 
# Ds12H3 files:
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H3out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/forward_Ds12H3.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H3out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/reverse_Ds12H3.bedgraph &
samtools view -b /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/sort_Ds12H3out.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/nt1598/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/nt1598/data_bioinf/lab2/bowtie1020/Ds12H3.bedgraph 

