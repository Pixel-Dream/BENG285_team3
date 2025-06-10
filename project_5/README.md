# BENG 285 / BNFO 285 / ECE 204 - Team 3 Repository for Project 3

Welcome to the shared codebase for Team 3 in the BENG 285 / BNFO 285 / ECE 204 course.

This repository contains the code developed for Project 5.

## Team Members

- Aleysha Chen
- Banghua Xu
- Eric Xue
- Haowen Zhou
- Jiaming Weng
- Peiyuan Han

## How to run the code
Please firstly go to the project_5 directory and run the following command:

```bash
cd project_5
gdown --folder https://drive.google.com/drive/folders/18shpLHrWOG4VzCYHpq9j9wodZcyc5yHk?usp=sharing
```

This will download the data to the `data` directory.

The mutation data is processed in the `mutation_feature_selection.Rmd` file. The further data preprocessing code is in the `dev_codes/data_pp.ipynb` file. We've provided the processed dataset from previous gdown session so no need to rerun it. The data is from various previous projects.

For SVM model training, cross-validation and hyperparam tuning, please refer to `linearSVM_ElasticNet.ipynb` file.

For Neural Network model training, cross-validation and hyperparam tuning, please run the following command:
```bash
# For continuous model
python train_nn_cv.py --model_type continuous --param_grid_json param_grid.json --output_csv continuous_nn_ablation_results.csv
# For Attention model
python train_nn_cv.py --model_type attention --param_grid_json param_grid.json --output_csv attention_nn_ablation_results.csv
```

For XGBoost model training, cross-validation and hyperparam tuning, please refer to `xgboost.ipynb` file.

For LLaMA-3.2-3B model, please firstly run the `llama3_data_cleaning.py` to preprocess the data. Then, please refer to `llama3_train.py` to train the model. We do not suggest running the `llama3_train.py` file because it takes too much resources to train the model.

For the ablation tests on SVM model, please run the following command:
```bash
python train_svm.py
```