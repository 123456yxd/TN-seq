#!/bin/bash
#SBATCH -J WT
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH -o WT.out


mkdir -p 01_clean 02_bam 03_counts
#index= ~/Reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.toplevel.fa.gz
gtf= ~/Reference/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.46.gtf

for sample in WT1 WT2 WT3 TN1 TN2 TN3; do
    # 1. 质控
    fastp -i raw/${sample}_1.fq.gz -I raw/${sample}_2.fq.gz \
          -o 01_clean/${sample}_R1.fq.gz -O 01_clean/${sample}_R2.fq.gz \
          -g -q 5 -u 50 -n 15 -l 150 \
          --overlap_diff_limit 1 --overlap_diff_percent_limit 10 \
          --detect_adapter_for_pe \
          -j 01_clean/${sample}.json -h 01_clean/${sample}.html -w 8

    # 2. 比对 + 排序
    bwa mem  ~/Reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.toplevel.fa.gz 01_clean/${sample}_R1.fq.gz 01_clean/${sample}_R2.fq.gz > ${sample}.sam
    samtools sort -@ 8 -o 02_bam/${sample}.sorted.bam ${sample}.sam
    samtools index 02_bam/${sample}.sorted.bam
done

# 3. 一次性定量
featureCounts -T 8 -p -t exon -g gene_id -a $gtf -o 03_counts/all_counts.txt 02_bam/*.sorted.bam
