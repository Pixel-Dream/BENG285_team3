# BENG 285 / BNFO 285 / ECE 204 - Team 3 Repository for Project 3

Welcome to the shared codebase for Team 3 in the BENG 285 / BNFO 285 / ECE 204 course.

This repository contains the code developed for Project 3.

## Team Members

- Aleysha Chen
- Banghua Xu
- Eric Xue
- Haowen Zhou
- Jiaming Weng
- Peiyuan Han

## How to run the code
Please firstly go to the project_4 directory and run the following command:

```bash
cd project_4
gdown --folder https://drive.google.com/drive/folders/1HhJ02YQFGrVWEo23duM-T_HBGiSMz-oL?usp=drive_link
```

This will download the data to the `data` directory.

Then, please run the following command to preprocess the metadata:

```bash
Rscript preprocessing.R
```

Then, for the DESeq2 analysis, please run the following command:

```bash
python deseq2.py
```

The results for DESeq2 will be saved to the `deseq2_corrected_age_gender_race` directory. And are visualized through the `plotVolcano.R` script.

For the driver gene analysis, please refer to the `driver_gene_mutation_MI.ipynb` file.

For the survival analysis, we have two versions of the code, including `proj4.R` and `Penalized Cox survival.ipynb`.
