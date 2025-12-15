# Human-genomics
 - This project focuses on variant analysis of human chromosome 7 specifically targeting the CFTR gene associated with Cystic Fibrosis.
 - The workflow uses raw FASTQ files, aligns them to the hg38 chromosome 7 reference and performs variant calling using GATK HaplotypeCaller, and annotates variants related to CFTR.
 ## Tools:
  - FastQC
  - BWA-MEM
  - SAMTools
  - GATK
    
## Documentation
- The complete workflow documentation is available in the Documentation folder
- [01. Introduction of this project](Documentation/01_Introduction.md)
- [02. Setup and Reference Genome Preparation](Documentation/02_Setup_and_Reference.md)
- [03. FastQC files](Documentation/03_FastQC_File.md)
- [04. Read Alignment](Documentation/04_Read_Alignment.md)
- [05. Marking Duplicates](Documentation/05_Marking_Duplicates.md)
