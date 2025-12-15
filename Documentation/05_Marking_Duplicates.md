Marking Duplicate Reads (GATK MarkDuplicates)

During PCR amplification and sequencing, the same DNA fragment can be sequenced multiple times, producing duplicate reads. These duplicates can bias variant allele frequencies and affect downstream analyses.

GATK MarkDuplicates identifies such duplicate reads and marks them in the BAM file, allowing variant callers to ignore them during analysis.
