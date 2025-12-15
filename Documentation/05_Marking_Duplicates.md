# IV.Marking Duplicate Reads (GATK MarkDuplicates)
During PCR amplification and sequencing the same DNA fragment can be sequenced multiple times producing duplicate reads. 

These duplicates can bias variant allele frequencies and affect downstream analyses.

GATK MarkDuplicates identifies such duplicate reads and marks them in the BAM file, allowing variant callers to ignore them during analysis.
  
---
## IV.1 - Why Read Groups (RG) Are Required

 - *GATK tools such as MarkDuplicates, BaseRecalibrator (BQSR), and HaplotypeCaller require Read Group (RG) information in the BAM file*.

 - Read Groups allow GATK to:

     Identify which reads belong to which sample
 
     Distinguish data from different lanes, libraries, or flow cells
 
     Track sequencing metadata such as:Platform (ILLUMINA), Library name, Sample ID

 ## IV.2 - Adding Read Groups to the Sorted BAM File

 - Read groups are added using *GATK Add or Replace ReadGroups*.

**Command**
```bash
  sudo docker run \ -v /workspaces/Human-genomics:/data \ broadinstitute/gatk:latest gatk AddOrReplaceReadGroups \  -I /data/Output_Files/Alignment/aln.sorted.bam \
  -O /data/Output_Files/Alignment/aln.RG.bam \ -RGID 1 \ -RGLB lib1 \ -RGPL ILLUMINA \ -RGPU unit1 \ -RGSM SRR30615628
```
 - This command generates output file of *aln.RG.bam* (Sorted BAM file with read group information added)

## IV.3 - Index the Read Group BAM File

 - The BAM file must be indexed before further processing.

**Command**
```bash
samtools index Output_Files/Alignment/aln.RG.bam
```
- This command generates output file of *aln.RG.bam.bai*

## IV.4 - Mark Duplicate Reads

 - Once read groups are present duplicate reads can be identified and marked

**Command**
```bash
sudo docker run \
  -v /workspaces/Human-genomics:/data \  broadinstitute/gatk:latest gatk MarkDuplicates \  -I /data/Output_Files/Alignment/aln.RG.bam \
  -O /data/Output_Files/Alignment/dedup.bam \  -M /data/Output_Files/Alignment/dup_metrics.txt
```
- This command generates output files of *dedup.bam *(BAM file with duplicate reads flagged), *dup_metrics.txt* (Summary statistics of detected duplicates

## IV.5 - Index the Deduplicated BAM File

 - The final BAM file must be indexed for downstream analysis.

**Command**
```bash
samtools index Output_Files/Alignment/dedup.bam
```
- This command generates the output file of *dedup.bam.bai*

## IV.6 - Output files 

  
