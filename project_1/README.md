# BENG 285 / BNFO 285 / ECE 204 - Team 3 Repository for Project 1

Welcome to the shared codebase for Team 3 in the BENG 285 / BNFO 285 / ECE 204 course.

This repository contains the code developed for Project 1.

## Team Members

- Aleysha Chen
- Banghua Xu
- Eric Xue
- Haowen Zhou
- Jiaming Weng
- Peiyuan Han

## How to run the code

1. Clone the repository

```bash
git clone https://github.com/your-repo/BENG285_team3.git
git checkout project_1
```

2. Install the Python dependencies

```bash
pip install -r requirements.txt
```

3. Download the data from [here](https://drive.google.com/drive/folders/1kANjAFNLurRn_vOBN3vw-OWZVtAVt2eF) and put it in the `data` folder. Or run the following command to download it automatically.

```bash
gdown --folder https://drive.google.com/drive/folders/1kANjAFNLurRn_vOBN3vw-OWZVtAVt2eF
```

4. Run the clustering script
The Jupyter Notebooks for linear dimensionality reduction (PCA and MDS) are under `clustering/linear/`. You can run them directly.
The codebased for nonlinear dimensionality reduction (t-SNE and UMAP) are under `clustering/nonlinear/`. You may run through this bash command. Parameters are set under `main.py`.

```bash
cd clustering/nonlinear
python main.py
```

5. Run the scoring script
The code for scoring is under `scoring_stats/`. You may run through this bash command.

```bash
cd scoring_stats
python scoring.py
```

For UMAP and t-SNE with best parameters, you can run the Jupyter Notebooks `scoring_stats/best_score_results.ipynb`.

6. Get UMAP and t-SNE results with the re-downloaded TCGA data
You can run the Jupyter Notebooks `re_downloaded_data_viz.ipynb` to get the results.


