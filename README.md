# NGS_analysis_pipeline
**The whole SARS-CoV-2 NGS data processing pipeline consists of 3 parts. You can use it separately or together.**
The following tools are required (current versions): bbtools, bwa, samtools, bcftools, mafft. 

Python libraries are required: Pandas, Numpy, Bio.



## 1. fq2snp.py 

Takes pairs of fastq files (pair-read), removes low-quality reads and some adapters.
Obtains bam-files without duplicates.
Obtains SNP-calling results.
Obtains mpileups (list of bases for all positions) as txt files.

### Usage:

```
python fq2snp.py <reference.fasta> <reads.fastq> -ot <mpileup_folder> -os <snp_folder>
```

**<reference.fasta>** - path to SARS-CoV-2 reference fasta file. 

**<reads.fastq>** - path to paired fastq files. Example of filenames: Cov-d186dl-290_S1_L001_R1_001.fastq.gz, Cov-d186dl-290_S1_L001_R2_001.fastq.gz  The suffix '_R1'/'_R2' in title is required.

**<mpileup_folder>** - path to folder for saving mpileup files. Default: mpileups

**<snp_folder>** - path to folder for saving vcf files after SNP-calling. Default: snps

## 2. m2f.py 

Processes mpileup files to fasta sequences (only single chromosome). 

Obtains the one fasta file for all sequences.

### Usage:

```
python m2f.py <mpileups> -o <fasta_folder>
```

**<mpileups>** - path to mpileup txt files

**<fasta_folder>** - path to folder for saving output fasta file. Default: generated_fastas

## 3. table_snp.py

For all  positions with SNPs obtains number of reads, reference base, alternative base, fractions of each base and presence in the subset of GISAID database.

### Usage:

```
python table_snp.py <reference.fasta> <snps.vcf> -o <output_folder> -db <database_subset_folder>
```

**<reference.fasta>** - path to SARS-CoV-2 reference fasta file. 

**<snps.vcf>** - path to vcf files with SNP-calling results.

**<output_folder>** - path to output folder. Default: csv_results

**<database_subset_folder>** - path to folder with GISAID subset. Default: fastas

## 4. aq2al.py

Uses all 3 mentioned above scripts together to obtain bam files, SNP-calling result file, consensus fasta file and multiple alignment between obtained sequences in fasta files.

### Usage:

```
python aq2al.py <reads.fastq> <reference.fasta>
```

**<reads.fastq>** - path to paired fastq files. Example of filenames: Cov-d186dl-290_S1_L001_R1_001.fastq.gz, Cov-d186dl-290_S1_L001_R2_001.fastq.gz  The suffix '_R1'/'_R2' in title is required.

**<reference.fasta>** - - path to SARS-CoV-2 reference fasta file. 

## Example:

```
python aq2al.py covid_files/* ref/sars-ref.fasta
```

