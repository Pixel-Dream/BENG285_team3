# BENG 285 / BNFO 285 / ECE 204 - Team 3 Repository for Project 2

Welcome to the shared codebase for Team 3 in the BENG 285 / BNFO 285 / ECE 204 course.

This repository contains the code developed for Project 2.

## Team Members

- Aleysha Chen
- Banghua Xu
- Eric Xue
- Haowen Zhou
- Jiaming Weng
- Peiyuan Han

## How to run the code

### Install python dependencies
```bash
pip install -r requirements.txt
```

### Data Preprocessing (`data_preprocessing.py`)
To run:
```bash
python data_preprocessing.py
```
Make sure you have the input file `TCGA.LUAD.fully_annotated.txt` in the `data` directory and the necessary libraries (`pandas`, `pybiomart`) installed. You will get a `results/result_with_gene_length.csv` file as output. 

### Liftover (`liftover.py`)
To run:
```bash
python liftover.py -i data/TCGA.LUAD.mutations.txt -c data/hg19ToHg38.over.chain.gz -r data/GRCh38_full_analysis_set_plus_decoy_hla.fa -o data/TCGA.LUAD.mutations.hg38.txt
```
Make sure you have the input file `TCGA.LUAD.mutations.txt`, chain file `hg19ToHg38.over.chain.gz`, and reference fasta file `GRCh38_full_analysis_set_plus_decoy_hla.fa` in the `data` directory and the necessary libraries (`pyliftover`, `pysam`) installed. You will get a `data/TCGA.LUAD.mutations.hg38.txt` file as output. 

You can download the chain file and reference fasta file through the following commands:

```bash
cd data
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

### Simplified dNdSpy (`dNdSpy`)
Reference genome and annotation file are required to run the `dNdSpy`. You can download the reference genome and annotation file through the following commands:
```bash
cd data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz
```
To run:
```bash
cd dNdSpy
python run_analysis.py
```
### Simplified dNdSpy with AlphaMissense (`am_score.ipynb` and `am_testing.ipynb`)
You need to download the required files with the following commands to run the `am_score.ipynb` and `am_testing.ipynb` files.
```bash
cd data
mkdir alphamissense
cd alphamissense
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_isoforms_hg38.tsv.gz
```
You may run the `am_score.ipynb` file to retrieve the AM scores for each variant on gene list from Simplified dNdSpy, then we calculate the expected and observed sum of AM scores for each gene. Finally run the `am_testing.ipynb` file to get the observed vs. expected sum_am_score/ds_count scores ratio for each gene.

### Permutation FDR (`permutation_pipeline.R`)
The permutation method is implemented in the `permutation_pipeline.R` file. And there is a function depo called `permutation_function.R` that contains the function used in the `permutation_pipeline.R` file.

To run:
```bash
Rscript permutation_pipeline.R path_to_annotation output_path 2000
```

### Benchmark (`benchmark.Rmd`)
Please download the following data (`IntOGen-DriverGenes_LUAD.tsv`, `dndscv_gene_level.csv`, `gene_results_with_am_significance.csv`) with the following commands:
```
bash
cd data
mkdir benchmark
gdown https://drive.google.com/file/d/1rbs4S81qbHGBntmRN_4dL-zY8omBQt_D/view?usp=drive_link
gdown https://drive.google.com/file/d/1aGiFgcbvYrARFfDqakZ1uj-XWkoCUKN6/view?usp=drive_link
gdown https://drive.google.com/file/d/1wUKPJ8sNZdwSmJf9T3O_gUs08RFzoX39/view?usp=drive_link
```
The `gene_results_with_am_significance.csv` file is the output of the `alphamissense` part. It would take a very long time to run the `alphamissense` part, so we provide the result file here.

This benchmark is mainly for plotting the results.

### Evaluation of statistical significance (`evaluation.ipynb`)
Please download the following data 
```bash
cd data
mkdir evaluation
gdown https://drive.google.com/file/d/1px_ibOZvY1gxGfKqwOIqcpoHX5uuT4x3/view?usp=drive_link
gdown https://drive.google.com/file/d/1PH--kwCnnYu2Lb9FfJCupbAvkAy11Onv/view?usp=drive_link
gdown https://drive.google.com/file/d/1EB-WQv3oje37Q6-GxMQ2mJsYH-HwewV_/view?usp=drive_link
```
The `evaluation.ipynb` file can be run to evaluate statistical significance of the results and used to directly visualize nonsynonymous mutation q-values vs. driver mutation posterior probabilities under mixture models. This file also contains code to plot a Venn diagram of the overlap between the driver genes identified by listed approaches.


