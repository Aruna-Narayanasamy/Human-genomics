#!/bin/bash

# ============================================================
# Complete GATK Variant Analysis Pipeline – Chromosome 7 (hg38)
# ============================================================
# Sample      : SRR30615628
# Reference   : chr7.fna
# Tools       : FastQC, BWA, SAMtools, GATK (Docker)
# ============================================================

set -euo pipefail

# -------------------------
# Create directories
# -------------------------
mkdir -p Reference_Genome
mkdir -p Raw_Data
mkdir -p Output_Files
mkdir -p Known_Sites
mkdir -p Output_Files/QC_Files
mkdir -p Output_Files/Alingnment
mkdir -p Output_Files/BQSR
mkdir -p Output_Files/Variant_Calling
mkdir -p Output_Files/Variant_Filtering

echo "Directory structure ready."

# ==========================================================
# STEP 1: Download Reference Genome and FASTQ Files
#===========================================================

echo "Downloding Data..."

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr7.fna.gz
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-run-15/SRR30615628/SRR30615628_1.fastq.gz
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-run-15/SRR30615628/SRR30615628_2.fastq.gz

# ============================================================
# STEP 1: Reference Genome Preparation
# ============================================================

echo "Unzipping reference genome..."

gunzip Reference_Genome/chr7.fna.gz

echo "Indexing reference genome..."

samtools faidx Reference_Genome/chr7.fna

echo "Creating sequence dictionary..."

sudo docker run \
-v /workspaces/Human-genomics:/data \
broadinstitute/gatk:latest gatk CreateSequenceDictionary \
-R /data/Reference_Genome/chr7.fna \
-O /data/Reference_Genome/chr7.dict
  

# ============================================================
# STEP 2: Quality Control (FastQC)
# ============================================================

echo "Running FastQC..."

fastqc \
 Raw_Data/SRR30615628_1.fastq.gz \ Raw_Data/SRR30615628_2.fastq.gz \ -o Output_Files/QC_Files/
 
# ============================================================
# STEP 3: Alignment using BWA-MEM
# ============================================================

echo "Indexing reference for BWA..."
bwa index  Reference_Genome/chr7.fna

echo "Running BWA-MEM alignment..."

bwa mem Reference_Genome/chr7.fna \
    Raw_Data/SRR30615628_1.fastq.gz \
    Raw_Data/SRR30615628_2.fastq.gz \
    > Output_Files/Alignment/aln.sam

# ============================================================
# STEP 4: SAM → BAM → Sort → Index
# ============================================================

echo "Converting SAM to sorted BAM..."
samtools view -S -b Output_Files/Alignment/aln.sam \ > Output_Files/Alignment/aln.bam
samtools sort Output_Files/Alignment/aln.bam \ -o Output_Files/Alignment/aln.sorted.bam
samtools index Output_Files/Alignment/aln.sorted.bam

rm ${ALIGN_DIR}/aln.sam ${ALIGN_DIR}/aln.bam

# ============================================================
# STEP 5: Add Read Groups
# ============================================================

echo "Adding read groups..."

sudo docker run \
 -v /workspaces/Human-genomics:/data \ 
 broadinstitute/gatk:latest gatk AddOrReplaceReadGroups \ 
  -I /data/Output_Files/Alignment/aln.sorted.bam \
  -O /data/Output_Files/Alignment/aln.RG.bam \ 
  -RGID 1 \ 
  -RGLB lib1 \ 
  -RGPL ILLUMINA \ 
  -RGPU unit1 \
  -RGSM SRR30615628

echo "Index the read group..."

samtools index Output_Files/Alignment/aln.RG.bam

# ============================================================
# STEP 6: Mark Duplicates
# ============================================================

echo "Marking duplicate reads..."

sudo docker run \
  -v /workspaces/Human-genomics:/data \  
  broadinstitute/gatk:latest gatk MarkDuplicates \  
  -I /data/Output_Files/Alignment/aln.RG.bam \
  -O /data/Output_Files/Alignment/dedup.bam \ 
  -M /data/Output_Files/Alignment/dup_metrics.txt

echo " Index the duplicate bam file..."
samtools index Output_Files/Alignment/dedup.bam

#==========================================================
# STEP 7: Known sites preparation 
#==========================================================

echo "Downloading known variants..."

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

echo "Restrict Known Sites to Chromosome 7..."

bcftools view -r 7 \
  Known_Sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \ -Oz -o Known_Sites/Mills_chr7.vcf.gz

echo "Index the Chromosome 7 Known Sites VCF..."

sudo docker run \
  -v $(pwd):/data \ broadinstitute/gatk:latest gatk IndexFeatureFile \ -I /data/Known_Sites/Mills_chr7.vcf.gz

echo "Convert .csi Index to .tbi..."

tabix -f -p vcf Known_Sites/Mills_chr7.vcf.gz

# ============================================================
# STEP 8: Base Quality Score Recalibration (BQSR)
# ============================================================

echo "Running BaseRecalibrator..."

sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk BaseRecalibrator \
  -R /data/Reference_Genome/chr7.fna \
  -I /data/Output_Files/Alignment/dedup.bam \
  --known-sites /data/Known_Sites/Mills_chr7.vcf.gz \
  -O /data/Output_Files/BQSR/recal_data.table

echo "Applying BQSR..."

sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk ApplyBQSR \
  -R /data/Reference_Genome/chr7.fna \
  -I /data/Output_Files/Alignment/dedup.bam \
  --bqsr-recal-file /data/Output_Files/BQSR/recal_data.table \
  -O /data/Output_Files/BQSR/aln.recal.bam

echo "Indexing recalibrated bam file..."

samtools index Output_Files/BQSR/aln.recal.bam

# ============================================================
# STEP 9: Variant Calling (HaplotypeCaller)
# ============================================================

echo "Running HaplotypeCaller..."

sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk HaplotypeCaller \
  -R /data/Reference_Genome/chr7.fna \
  -I /data/Output_Files/BQSR/aln.recal.bam \
  -O /data/Output_Files/Variant_Calling/chr7_variants_bqsr.vcf.gz

# ============================================================
# STEP 10: Variant Separation
# ============================================================

echo "Separating SNPs and INDELs..."

sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk SelectVariants \
  -R /data/Reference_Genome/chr7.fna \
  -V /data/Output_Files/Variant_Calling/chr7_variants_bqsr.vcf.gz \
  --select-type-to-include SNP \
  -O /data/Output_Files/Variant_Filtering/chr7_snps.vcf.gz

sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk SelectVariants \
  -R /data/Reference_Genome/chr7.fna \
  -V /data/Output_Files/Variant_Calling/chr7_variants_bqsr.vcf.gz \
  --select-type-to-include INDEL \
  -O /data/Output_Files/Variant_Filtering/chr7_indels.vcf.gz

# ============================================================
# STEP 10: Variant Filtration
# ============================================================

echo "Filtering SNPs..."

sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk VariantFiltration \
  -R /data/Reference_Genome/chr7.fna \
  -V /data/Output_Files/Variant_Filtering/chr7_snps.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  -O /data/Output_Files/Variant_Filtering/chr7_snps_filtered.vcf.gz

echo "Filtering INDELs..."

sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk VariantFiltration \
  -R /data/Reference_Genome/chr7.fna \
  -V /data/Output_Files/Variant_Filtering/chr7_indels.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "SOR > 10.0" --filter-name "SOR10" \
  -O /data/Output_Files/Variant_Filtering/chr7_indels_filtered.vcf.gz

 echo "Retain Only High-Confidence (PASS) Variants..."
 echo "SNPs..."

 sudo docker run \
  -v $(pwd):/data \ 
  broadinstitute/gatk:latest gatk SelectVariants \ 
  -V /data/Output_Files/Variant_Filtering/chr7_snps_filtered.vcf.gz \
  --exclude-filtered \ 
  -O /data/Output_Files/Variant_Filtering/chr7_snps_pass.vcf.gz

  echo "INDELs..."

  sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk SelectVariants \ -V /data/Output_Files/Variant_Filtering/chr7_indels_filtered.vcf.gz \
  --exclude-filtered \ 
  -O /data/Output_Files/Variant_Filtering/chr7_indels_pass.vcf.gz

  echo "Merge PASS SNPs and INDELs..."

  sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk MergeVcfs \
  -I /data/Output_Files/Variant_Filtering/chr7_snps_pass.vcf.gz \
  -I /data/Output_Files/Variant_Filtering/chr7_indels_pass.vcf.gz \
  -O /data/Output_Files/Variant_Filtering/chr7_variants_filtered_final.vcf.gz

echo "Pipeline completed successfully."
