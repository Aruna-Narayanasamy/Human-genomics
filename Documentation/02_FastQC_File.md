# II. FASTQ Quality Check (FastQC)

Quality control is performed on the raw paired-end FASTQ files using **FastQC**.  
This step helps assess:

- Read quality  
- Adapter contamination  
- GC content  
- Overrepresented sequences  
- Any potential sequencing issues  

---

## II.1 – Install FastQC

 FastQC installed using `apt`:

```bash
sudo apt-get install -y fastqc
```

## II.2 – Run FastQC on FASTQ Files

 -FastQC is executed on both paired-end FASTQ files.

**1.Create output folder**
 An output folder is created to store all resulting files.
 Inside the main Output_Files directory, a subfolder named QC_Files is created specifically for storing the FastQC reports.

-**Command**:

```bash
 mkdir Output_Files
 mkdir -p Output_Files/QC_Files
```
**2.Run FastQC**

-**Command**

```bash
   fastqc Raw_Data/SRR30615628_1.fastq.gz \ Raw_Data/SRR30615628_2.fastq.gz \ -o Output_Files/QC_Files/
```
 FastQC produce two files of each FASTQ file (.html, .zip).

## II.3 - Files generated in this step

Output_Files/QC_Files/
- SRR30615628_1_fastqc.html
- SRR30615628_1_fastqc.zip
- SRR30615628_2_fastqc.html
- SRR30615628_2_fastqc.zip
