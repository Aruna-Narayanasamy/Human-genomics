# VII. Variant Calling Using GATK HaplotypeCaller – Chromosome 7
This document describes *variant calling* using *GATK HaplotypeCaller* on a *BQSR-recalibrated BAM file* for *human chromosome 7 (hg38)*.

Variant calling identifies genomic variants such as *Single Nucleotide Polymorphisms (SNPs)* and *insertions/deletions (INDELs)* by comparing aligned
sequencing reads against a reference genome.

In GATK-based pipelines *HaplotypeCaller* performs local de-novo assembly of haplotypes and calls variants with high accuracy. It is executed after all
preprocessing steps including *duplicate marking and BQSR*.

---

## Input Files needed

| File | Description |
|----|------------|
| `Reference_Genome/chr7.fna` | Reference genome (chromosome 7) |
| `Output_Files/BQSR/aln.recal.bam` | BQSR-recalibrated BAM file |
| `Output_Files/BQSR/aln.recal.bam.bai` | BAM index file |

---

## VII.1 - Create Variant Calling Output Directory

- Creates a dedicated directory to store variant calling results.

  **Command**
```bash
mkdir -p Output_Files/Variant_Calling
```
## VII.2 - Run GATK HaplotypeCaller

- Variant calling was performed using the recalibrated BAM file as input.

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \ broadinstitute/gatk:latest gatk HaplotypeCaller \  -R /data/Reference_Genome/chr7.fna \
  -I /data/Output_Files/BQSR/aln.recal.bam \  -O /data/Output_Files/Variant_Calling/chr7_variants_bqsr.vcf.gz
```
- This command generates the file of chr7_variants_bqsr.vcf.gz, chr7_variants_bqsr.vcf.gz.tbi

## VII.3 - Files generated in this step
 - Output_Files/Variant_Calling
     - chr7_variants_bqsr.vcf.gz
     - chr7_variants_bqsr.vcf.gz.tbi
