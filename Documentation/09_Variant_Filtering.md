# VIII. Variant Filtering – Chromosome 7

This document describes the *variant filtering* steps performed after variant calling using *GATK HaplotypeCaller*.

Filtering was applied to retain *high-confidence SNPs and INDELs* for *human chromosome 7 (hg38)*.

Variant calling produces both *true variants* and *potential false positives*.

Variant filtering applies quality-based criteria to distinguish *high-confidence variants* from low-quality calls.

Because this analysis involves a *single-sample, chromosome-level dataset*, *hard filtering* was used instead of VQSR.

## Input files needed for variant calling

- Reference_Genome/chr7.fna
- Variant_Calling/chr7_variants_bqsr.vcf.gz

---

## VIII.1 - Create Variant Filtering Output Directory

**Command**
```bash
mkdir -p Output_Files/Variant_Filtering
```
## VIII.2 - Separate SNPs and INDELs
- SNPs and INDELs are separated because they have different error profiles and require different filtering thresholds.

### VIII.2.1 - Extract SNPs

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \ broadinstitute/gatk:latest gatk SelectVariants \  -R /data/Reference_Genome/chr7.fna \
  -V /data/Output_Files/Variant_Calling/chr7_variants_bqsr.vcf.gz \ --select-type-to-include SNP \
  -O /data/Output_Files/Variant_Filtering/chr7_snps.vcf.gz
```
- This command generates the file of *chr7_snps.vcf.gz, chr7_snps.vcf.gz.tbi*

### VIII.2.2 - Extract INDELs

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \  broadinstitute/gatk:latest gatk SelectVariants \  -R /data/Reference_Genome/chr7.fna \
  -V /data/Output_Files/Variant_Calling/chr7_variants_bqsr.vcf.gz \  --select-type-to-include INDEL \
  -O /data/Output_Files/Variant_Filtering/chr7_indels.vcf.gz
```
- This command generates the file of *chr7_indels.vcf.gz, chr7_indels.vcf.gz.tbi*

## VIII.3 - Apply Hard Filters to SNPs

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \  broadinstitute/gatk:latest gatk VariantFiltration \-R /data/Reference_Genome/chr7.fna \
-V/data/Output_Files/Variant_Filtering/chr7_snps.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  -O /data/Output_Files/Variant_Filtering/chr7_snps_filtered.vcf.gz
```
- This command generates the file of *chr7_snps_filtered.vcf.gz, chr7_snps_filtered.vcf.gz.tbi*
- Filter Criteria Used:
  - QD (Quality by Depth)
  - FS (Fisher Strand Bias)
  - MQ (Mapping Quality)
  - SOR (Strand Odds Ratio

## VIII.4 - Apply Hard Filters to INDELs

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \  broadinstitute/gatk:latest gatk VariantFiltration \  -R /data/Reference_Genome/chr7.fna \
-V/data/Output_Files/Variant_Filtering/chr7_indels.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "SOR > 10.0" --filter-name "SOR10" \
  -O /data/Output_Files/Variant_Filtering/chr7_indels_filtered.vcf.gz
```
- This command generates the file of *chr7_indels_filtered.vcf.gz, chr7_indels_filtered.vcf.gz.tbi*

## VIII.5 - Retain Only High-Confidence (PASS) Variants

- Only variants with FILTER=PASS are retained.
### VIII.5.1 - PASS SNPs

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \ broadinstitute/gatk:latest gatk SelectVariants \ -V /data/Output_Files/Variant_Filtering/chr7_snps_filtered.vcf.gz \
  --exclude-filtered \ -O /data/Output_Files/Variant_Filtering/chr7_snps_pass.vcf.gz
```
- This command generates the file of *chr7_snps_pass.vcf.gz, chr7_snps_pass.vcf.gz.tbi*

### VIII.5.2 - PASS INDELs

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk SelectVariants \ -V /data/Output_Files/Variant_Filtering/chr7_indels_filtered.vcf.gz \
  --exclude-filtered \ -O /data/Output_Files/Variant_Filtering/chr7_indels_pass.vcf.gz
```
- This command generates the file of *chr7_indels_pass.vcf.gz, chr7_indels_pass.vcf.gz.tbi*

## VIII.6 -  Merge PASS SNPs and INDELs

**Command**
```bash
sudo docker run \
  -v $(pwd):/data \
  broadinstitute/gatk:latest gatk MergeVcfs \ -I /data/Output_Files/Variant_Filtering/chr7_snps_pass.vcf.gz \
  -I /data/Output_Files/Variant_Filtering/chr7_indels_pass.vcf.gz \ -O /data/Output_Files/Variant_Filtering/chr7_variants_filtered_final.vcf.gz
```
- This command generates the file of *chr7_variants_filtered_final.vcf.gz, chr7_variants_filtered_final.vcf.gz.tbi*

## VIII.7 - Files generated in this step
 - Output_Files/Variant_Filtering
   - chr7_snps.vcf.gz
   - chr7_snps.vcf.gz.tbi
   - chr7_indels.vcf.gz
   - chr7_indels.vcf.gz.tbi
   - chr7_snps_filtered.vcf.gz
   - chr7_snps_filtered.vcf.gz.tbi
   - chr7_indels_filtered.vcf.gz
   - chr7_indels_filtered.vcf.gz.tbi
   - chr7_snps_pass.vcf.gz
   - chr7_snps_pass.vcf.gz.tbi
   - chr7_indels_pass.vcf.gz
   - chr7_indels_pass.vcf.gz.tbi
   -  chr7_variants_filtered_final.vcf.gz
   -  chr7_variants_filtered_final.vcf.gz.tbi
