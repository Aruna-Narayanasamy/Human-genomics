# VI. Base Quality Score Recalibration (BQSR)

BQSR analyzes patterns of covariation in the sequence dataset to correct systematic errors made by the sequencing machine.

---
## VI.1 - Create BQSR Output Directory

**Command**
```bash
mkdir -p Output_Files/BQSR
```
## VI.2 - Generate Recalibration Table (BaseRecalibrator)

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \ broadinstitute/gatk:latest gatk BaseRecalibrator \  -R /data/Reference_Genome/chr7.fna \  -I /data/Output_Files/Alignment/aln.RG.bam \
  --known-sites /data/Known_Sites/Mills_chr7.vcf.gz \ -O /data/Output_Files/BQSR/recal_data.table
```
- This command generates the file of *recal_data.table*

## VI.3 - Apply BQSR to BAM File

**Command**
```bash
-v $(pwd):/data \  broadinstitute/gatk:latest gatk ApplyBQSR \  -R /data/Reference_Genome/chr7.fna \  -I /data/Output_Files/Alignment/aln.RG.bam \
  --bqsr-recal-file /data/Output_Files/BQSR/recal_data.table \ -O /data/Output_Files/BQSR/aln.recal.bam
```
- This command generates the file of *aln.recal.bam*

## VI.4 - Index the Recalibrated BAM File

**Command**
```bash
samtools index Output_Files/BQSR/aln.recal.bam
```
- This command generates the file of aln.recal.bam.bai

## VI.5 - Files generated in this step
 - Output_Files/BQSR
    - recal_data.table
    - aln.recal.bam
    - aln.recal.bam.bai
   
