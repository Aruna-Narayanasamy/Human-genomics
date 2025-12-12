# III. Read Alignment Using BWA-MEM

Paired-end sequencing reads are aligned to the reference genome using **BWA-MEM**, a widely used algorithm for mapping short sequencing reads.
This step produces SAM and BAM files that are essential for downstream processes such as duplicate marking, base recalibration, and variant calling.

---

## III.1 Install BWA

Install BWA using apt:

```bash
sudo apt-get install -y bwa
```

## III.2 Create Output Directory for Alignment

All alignment-related files will be stored in Alignment folder.

**Command**:
```bash
mkdir -p Output_Files/Alignment
```

## III.3 Index the Reference Genome for BWA

BWA requires its own index files (different from .fai and .dict).

**Command**
```bash
bwa index Reference_Genome/chr7.fna
```
This command produces **chr7.fna.amb, chr7.fna.ann, chr7.fna.bwt, chr7.fna.pac, chr7.fna.sa**

## III.4 Run BWA-MEM Alignment

For this we need FASTQ files and reference genome fna file.

**Command**
```bash
bwa mem Reference_Genome/chr7.fna \
    Raw_Data/SRR30615628_1.fastq.gz \
    Raw_Data/SRR30615628_2.fastq.gz \
    > Output_Files/Alignment/aln.sam
```
This command produce **.sam file**

## III.5 Convert SAM to BAM

SAM files are large and uncompressed. Need to convert BAM format

**Command**
```bash
samtools view -S -b Output_Files/Alignment/aln.sam \ > Output_Files/Alignment/aln.bam
```
this cmmand generates **aln.bam file**

## III.6 Sort BAM File

Sorting is required for downstream tools such as GATK and samtools index.

**Command**:
```bash
samtools sort Output_Files/Alignment/aln.bam \ -o Output_Files/Alignment/aln.sorted.bam
```
This command generates **aln.sorted.bam**

## III.7 Index the Sorted BAM File

Create a BAM index (.bai) for random access during variant calling

**Command**
```bash
samtools index Output_Files/Alignment/aln.sorted.bam
```

This command generates **aln.sorted.bam.bai**

## III.8 - Files generated in this step
|Output_Files/Alignment| Reference_Genome|
|----------------------|-----------------|
|aln.sam |chr7.fna.|
|aln.bam |chr7.fna.ann|
|aln.sorted.bam|chr7.fna.bwt|
|aln.sorted.bam.bai|chr7.fna.pac|
|       |chr7.fna.sa|
