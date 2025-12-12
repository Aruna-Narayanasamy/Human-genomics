# Introduction

## **Human Genomics Project**

Human genomics is the study of the structure, function, evolution, and mapping of the human genome.  
By analyzing DNA sequences, researchers can understand genetic variation, identify disease-causing mutations, and develop targeted therapies.

This project demonstrates a complete workflow for processing human sequencing data, focusing on variant detection within a medically important genomic region.

---

## **What is Cystic Fibrosis?**

Cystic Fibrosis (CF) is a **genetic disorder** caused by mutations in the **CFTR gene**  
(*Cystic Fibrosis Transmembrane Conductance Regulator*).

- It affects **respiratory**, **digestive**, and **reproductive** systems.  
- The disease results from **defective chloride ion transport**, causing thick mucus buildup.  
- More than **2,000 CFTR variants** are known, but only some cause disease.  
- It is one of the most common life-shortening genetic disorders.

Since CF is a monogenic (single-gene) disorder, it is ideal for learning variant analysis.

## **Why Chromosome 7 Was Chosen**

Cystic Fibrosis is caused by variants in the **CFTR gene** which is located on:

- **Chromosome 7q31.2**
- Genomic coordinates (hg38): *chr7:117,120,017 – 117,308,718*

By working specifically with **chromosome 7**, this project becomes faster, more focused.

### Benefits of using only chr7:
- Reduced computational load  
- Faster alignment and variant calling  
- Easier to interpret CFTR-specific variants  
   
---

## **Project Objective**

This project aims to:

- Download real human sequencing data (FASTQ files)
- Prepare the reference genome (chr7)
- Run quality control on raw sequences
- Align reads using BWA-MEM
- Process alignment (sort, index, mark duplicates)
- Perform variant calling (HaplotypeCaller)
- Annotate variants (e.g., ClinVar, gene annotations)
- Identify potential CFTR gene variants

---

## **Workflow Overview**

The complete pipeline includes the following steps:

1. **Setup & Reference Genome Preparation**
   - Download chr7 reference
   - Index reference (FASTA index, sequence dictionary)

2. **Quality Control (FastQC)**
   - Assess sequencing quality and identify issues

3. **Read Alignment (BWA-MEM)**
   - Map FASTQ reads to the chr7 reference genome

4. **Post-Alignment Processing**
   - Convert SAM → BAM
   - Sort and index BAM
   - Add read groups
   - Mark duplicate reads (GATK MarkDuplicates)

5. **Variant Calling (GATK HaplotypeCaller)**
   - Generate VCF of SNPs and Indels

6. **Variant Filtering**
   - Apply quality filters to low-confidence variants

7. **Variant Annotation**
   - Interpret variants using tools like ANNOVAR or VEP
   - Identify CFTR-related mutations

---

## **Tools and Technologies Used**

| Tool / Software | Purpose |
|-----------------|---------|
| **FastQC** | FASTQ quality assessment |
| **BWA-MEM** | Read alignment |
| **samtools** | File conversion, sorting, indexing |
| **GATK** | MarkDuplicates, BQSR, variant calling |
| **Docker** | Running GATK efficiently |


---

## **Summary**

This project provides a practical, step-by-step human genomics workflow focusing on **chromosome 7** and the **CFTR gene**, which causes Cystic Fibrosis.  



