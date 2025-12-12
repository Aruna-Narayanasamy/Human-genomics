# I.Setup and Reference Genome preparation

  This section explains how to download the reference genome and raw FASTQ files, create project directories and prepare the reference genome for downstream analysis.

  ---
 
 ## I.1 - Download Reference Genome and FASTQ Files
 
  Raw sequencing data and the chromosome 7 reference genome are required to begin any alignment or variant-calling workflow.
  FASTQ Files (SRR30615628) contain paired-end sequencing reads from the SRA database.

 - **Command for Reference Genome**:
   
    ```bash
       wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr7.fna.gz
     ```
    
 - **Command for FASTQ**:
   
      ```bash
    
      wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-run-15/SRR30615628/SRR30615628_1.fastq.gz
      wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-run-15/SRR30615628/SRR30615628_2.fastq.gz
      ```

 ## I.2 - Create Project Directories
 
  I created folders to organize reference files and raw sequence data.
 
   - **Command**:
     
      ```bash
      mkdir Reference_Genome
      mkdir Raw_Data
      ```
   -  Reference_Genome/ - Stores chr7.fna, chr7.fna.fai, chr7.dict
   - Raw_Data/ - Stores FASTQ files (SRR30615628_1 & SRR30615628_2)

 ## I.3 - Move Input Files
   - **Command** :

   ```bash
     mv chr7.fna.gz Reference_Genome/
     mv SRR30615628_1.fastq.gz Raw_Data/
     mv SRR30615628_2.fastq.gz Raw_Data/
   ```

 ## I.4 - Reference Genome preparation
 
  Before alignment the reference genome must be prepared by unzipping, indexing, and creating a GATK sequence dictionary.
  
   **1.Unzip the Genome**
  
   -**Command**:

    ```bash
    gunzip Reference_Genome/chr7.fna.gz
    ```
  **2.Generate FASTA Index**

  The .fai index file enables fast random access to the reference sequence, required by tools like samtools, BWA, and GATK.
  
   -**Command**:

   ```bash
  samtools faidx Reference_Genome/chr7.fna
 ```
  **3.Create Sequence Dictionary (GATK + Docker)**

 GATK requires a dictionary file describing sequence names and lengths.
 I installed GATK using Docker, the command must run inside a container.

   -**Command**:

   ```bash
  sudo docker run \
  -v /workspaces/Human-genomics:/data \
  broadinstitute/gatk:latest gatk CreateSequenceDictionary \
  -R /data/Reference_Genome/chr7.fna \
  -O /data/Reference_Genome/chr7.dict
```

## I.5- Files generated in this step
| Folder            | Files                          |
|-------------------|--------------------------------|
| Reference_Genome  | `chr7.fna` <br> `chr7.fna.fai` <br> `chr7.dict` |
| Raw_Data          | `SRR30615628_1.fastq.gz` <br> `SRR30615628_2.fastq.gz` |


    



 
