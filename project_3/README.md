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

### Install python dependencies
```bash
pip install -r requirements.txt
```

### Mutational Matrix Generator (`matrix_generator/generate_sbs_matrix.py`)
To run:
```bash
bash matrix_generator/cmd.sh
```
Make sure you have the input file `TCGA.LUAD.mutations.txt` in the `data` directory and the necessary libraries (`pandas`, `pybiomart`) installed. You will get a `output/sbs_96_matrix.tsv` file as output. 

### Mutational Signature Extractor (`signature_extractor/extract_signatures.py`)
To run:
```bash
python signature_extractor/extract_signatures.py
```

### 