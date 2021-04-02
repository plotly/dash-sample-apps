
# Data Usage

All data is available in [GEO accession GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450).
We used mus musculus RNA-seq data from basal and luminal cell types (2 day lactate),
resulting in 4 samples (2 replicates per condition).

Files and processing was based off of an [RNA-seq tutorial from WEHI bioinformatics](https://bioinformatics-core-shared-training.github.io/RNAseq-R/). In all samples, we only evaluate reads from chr1.

# Data generation

4 samples was downloaded using SRA (SRP045534) using 'prefetch' in sra-toolkit.
Fastq files were extracted using 'fastq-dump' in sra-toolkit.

Fastq files were aligned to the mm10 assembly using Rsubread in R version 4.0.3.
Reads were quantified using `featureCounts()` in Rsubread.
Limma and EdgeR were used to compute differential expression between genes.

Bam files were then filtered using the following criteria:
- bam files were sorted using `samtools sort`
- bam files were filtered to only include regions overlapping genes (samtools view)
- regions were sampled to only include 20% of original reads (samtools view -s 0.2)
- final filtered bam files were indexed
- sometimes samtools filter generates invalid bam files. In this case, we run
		`samtools view` to convert to a sam file, then convert back to a bam file.

Entrez GeneIDs were translated to mgi gene symbols using biomaRt in R.
