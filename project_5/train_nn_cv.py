import argparse
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader
import json

from data_loader import load_and_preprocess_data
from model.models import ContinuousFeatureNN, AttentionFeatureNN
from training_utils import OmicsDataset, train_epoch, evaluate_model, tune_pytorch_hyperparameters

def construct_ablation_feature_subsets(adata, expression_type, mutation_type):
    """
    Constructs feature names for a specific ablation scenario.
    """
    all_rna_features = [f'gene_{g}' for g in adata.var_names]
    
    rna_features = []
    if expression_type:
        rna_features = all_rna_features

    dna_features = []
    if mutation_type:
        if mutation_type not in adata.uns:
            raise ValueError(f"Mutation key '{mutation_type}' not found in adata.uns")
        dna_features = [f'mut_{c}' for c in adata.uns[mutation_type].columns]

    return dna_features, rna_features

def get_feature_matrix_and_labels(adata, expression_type, mutation_type):
    """
    Constructs the full feature matrix and returns it with labels.
    This is a simplified version of logic from other scripts.
    """
    sample_ids = adata.obs_names
    
    expr_df = None
    if expression_type:
        if expression_type == 'normalized':
            data = adata.X
        elif expression_type == 'raw':
            data = adata.layers['counts']
        else:
            raise ValueError(f"Unknown expression type: {expression_type}")
        
        if hasattr(data, 'toarray'):
            data = data.toarray()
        
        expr_df = pd.DataFrame(data, index=sample_ids, columns=[f'gene_{g}' for g in adata.var_names])

    mut_df = None
    if mutation_type:
        if mutation_type not in adata.uns:
            raise ValueError(f"Mutation key '{mutation_type}' not found in adata.uns")
        mut_data = adata.uns[mutation_type].copy()
        mut_df = pd.DataFrame(mut_data, index=mut_data.index, columns=[f'mut_{c}' for c in mut_data.columns])

    if expr_df is not None and mut_df is not None:
        feature_matrix = expr_df.join(mut_df, how='left').fillna(0)
    elif expr_df is not None:
        feature_matrix = expr_df
    elif mut_df is not None:
        feature_matrix = mut_df.reindex(sample_ids).fillna(0)
    else:
        feature_matrix = pd.DataFrame(index=sample_ids)
        
    labels = adata.obs['survival_3yr_label'].values
    return feature_matrix, labels

def train_final_model(model, X_dna_train, X_rna_train, y_train, n_epochs, batch_size, learning_rate):
    """Train the model on the full training data."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    
    scaler_dna = None
    if X_dna_train.shape[1] > 0:
        scaler_dna = StandardScaler().fit(X_dna_train)
        X_dna_train = scaler_dna.transform(X_dna_train)

    scaler_rna = None
    if X_rna_train.shape[1] > 0:
        scaler_rna = StandardScaler().fit(X_rna_train)
        X_rna_train = scaler_rna.transform(X_rna_train)

    train_dataset = OmicsDataset(X_dna_train, X_rna_train, y_train)
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    
    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    
    for epoch in range(n_epochs):
        train_loss, train_acc = train_epoch(model, train_loader, criterion, optimizer, device)
        print(f"Epoch {epoch+1}/{n_epochs} - Train Loss: {train_loss:.4f}, Train Acc: {train_acc:.4f}")
        
    return model, scaler_dna, scaler_rna

def evaluate_final_model(model, X_dna_test, X_rna_test, y_test, scaler_dna, scaler_rna, batch_size):
    """Evaluate the final model on the hold-out test set."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    
    if scaler_dna:
        X_dna_test = scaler_dna.transform(X_dna_test)
    
    if scaler_rna:
        X_rna_test = scaler_rna.transform(X_rna_test)
    
    test_dataset = OmicsDataset(X_dna_test, X_rna_test, y_test)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    criterion = nn.CrossEntropyLoss()
    
    # Reusing evaluate_model from training_utils
    test_loss, test_auc = evaluate_model(model, test_loader, criterion, device)
    
    # We need other metrics too
    model.eval()
    all_labels = []
    all_preds = []
    with torch.no_grad():
        for batch in test_loader:
            dna = batch['dna'].to(device)
            rna = batch['rna'].to(device)
            labels = batch['label'].to(device)
            outputs = model(dna, rna)
            preds = torch.argmax(outputs, dim=1)
            all_labels.extend(labels.cpu().numpy())
            all_preds.extend(preds.cpu().numpy())

    return {
        'auroc': test_auc,
        'accuracy': accuracy_score(all_labels, all_preds),
        'balanced_accuracy': balanced_accuracy_score(all_labels, all_preds),
        'f1_score': f1_score(all_labels, all_preds, average='binary')
    }

def main(args):
    """Main script for hyperparameter tuning and ablation study."""
    # 1. Load data
    print("Loading data...")
    adata = load_and_preprocess_data(args.anndata_path)

    # 2. Outer train-test split (80% for CV, 20% for hold-out)
    indices = np.arange(adata.n_obs)
    train_indices, test_indices = train_test_split(
        indices,
        test_size=args.test_size,
        random_state=42,
        stratify=adata.obs['survival_3yr_label'].values
    )
    train_adata = adata[train_indices, :]
    test_adata = adata[test_indices, :]
    print(f"Full dataset size: {len(adata)}")
    print(f"Training/CV set size: {len(train_adata)}")
    print(f"Hold-out test set size: {len(test_adata)}")

    # 3. Model selection
    if args.model_type == 'continuous':
        model_class = ContinuousFeatureNN
    elif args.model_type == 'attention':
        model_class = AttentionFeatureNN
    else:
        raise ValueError(f"Unknown model type: {args.model_type}")

    # 4. Define hyperparameter grid
    with open(args.param_grid_json, 'r') as f:
        param_grid = json.load(f)

    # 5. Define ablation study scenarios
    ablation_scenarios = [
        {'name': 'RNA (normalized)', 'expression': 'normalized', 'mutation': None},
        {'name': 'RNA (raw)', 'expression': 'raw', 'mutation': None},
        {'name': 'Mutations (full)', 'expression': None, 'mutation': 'jiaming_data_full'},
        {'name': 'RNA (norm) + Mut (full)', 'expression': 'normalized', 'mutation': 'jiaming_data_full'},
        {'name': 'RNA (raw) + Mut (full)', 'expression': 'raw', 'mutation': 'jiaming_data_full'},
    ]
    
    # 6. Run ablation study
    print(f"\n--- Starting Ablation Study for {args.model_type} model ---")
    ablation_results = []

    for scenario in ablation_scenarios:
        print(f"\n--- Running scenario: {scenario['name']} ---")

        # Get feature names for this scenario
        dna_feature_names, rna_feature_names = construct_ablation_feature_subsets(
            train_adata, scenario['expression'], scenario['mutation']
        )
        
        feature_names = dna_feature_names + rna_feature_names
        if not feature_names:
            print(f"Skipping scenario '{scenario['name']}' due to having 0 features.")
            continue
            
        # Construct feature matrix from the training AnnData
        feature_matrix, labels = get_feature_matrix_and_labels(
            train_adata, scenario['expression'], scenario['mutation']
        )
        
        # Check for empty feature matrix
        if feature_matrix.shape[1] == 0:
            print(f"Skipping scenario '{scenario['name']}' due to having 0 features.")
            continue
            
        print(f"Feature matrix shape for tuning: {feature_matrix.shape}")

        static_model_kwargs = {
            'dna_feature_names': dna_feature_names,
            'rna_feature_names': rna_feature_names
        }

        # a. Hyperparameter search on the 80% data
        tuning_results = tune_pytorch_hyperparameters(
            feature_matrix=feature_matrix,
            labels=labels,
            model_class=model_class,
            param_grid=param_grid,
            static_model_kwargs=static_model_kwargs,
            k=args.k_folds,
            n_epochs=args.n_epochs,
            batch_size=args.batch_size
        )
        best_params = tuning_results['best_params']
        best_cv_score = tuning_results['best_score']
        all_tuning_results = pd.DataFrame(tuning_results['all_results'])
        
        print(f"Best parameters found: {best_params} with CV AUROC: {best_cv_score:.4f}")

        # Save detailed tuning results for this scenario
        if args.output_csv:
            scenario_name_safe = "".join(c for c in scenario['name'] if c.isalnum() or c in (' ', '_')).rstrip().replace(' ', '_')
            tuning_filename = args.output_csv.replace('.csv', f'_{scenario_name_safe}_tuning_results.csv')
            all_tuning_results.sort_values('avg_auc', ascending=False, inplace=True)
            all_tuning_results.to_csv(tuning_filename, index=False)
            print(f"Saved detailed tuning results to {tuning_filename}")

        # b. Train final model on the full 80% data with best params
        print("Training final model on full training set...")
        final_model_kwargs = {
            'n_dna_features': len(dna_feature_names),
            'n_rna_features': len(rna_feature_names),
            'hidden_dim': best_params['hidden_dim'],
            'dropout_rate': best_params['dropout_rate']
        }
        final_model = model_class(**final_model_kwargs)
        
        X_dna_train = feature_matrix[dna_feature_names].values if dna_feature_names else np.empty((len(feature_matrix), 0), dtype=np.float32)
        X_rna_train = feature_matrix[rna_feature_names].values if rna_feature_names else np.empty((len(feature_matrix), 0), dtype=np.float32)

        final_model, scaler_dna, scaler_rna = train_final_model(
            model=final_model,
            X_dna_train=X_dna_train,
            X_rna_train=X_rna_train,
            y_train=labels,
            n_epochs=args.n_epochs,
            batch_size=args.batch_size,
            learning_rate=best_params['learning_rate']
        )

        # c. Evaluate on the 20% hold-out test set
        print("Evaluating final model on hold-out test set...")
        test_feature_matrix, test_labels = get_feature_matrix_and_labels(
            test_adata, scenario['expression'], scenario['mutation']
        )
        
        X_dna_test = test_feature_matrix[dna_feature_names].values if dna_feature_names else np.empty((len(test_feature_matrix), 0), dtype=np.float32)
        X_rna_test = test_feature_matrix[rna_feature_names].values if rna_feature_names else np.empty((len(test_feature_matrix), 0), dtype=np.float32)

        test_metrics = evaluate_final_model(
            model=final_model,
            X_dna_test=X_dna_test,
            X_rna_test=X_rna_test,
            y_test=test_labels,
            scaler_dna=scaler_dna,
            scaler_rna=scaler_rna,
            batch_size=args.batch_size
        )
        
        print(f"Test Set Metrics: {test_metrics}")

        # Store results
        result_row = {
            'Scenario': scenario['name'],
            'Num DNA Features': len(dna_feature_names),
            'Num RNA Features': len(rna_feature_names),
            'best_cv_auroc': best_cv_score,
            **{f'test_{k}': v for k, v in test_metrics.items()},
            'best_params': str(best_params)
        }
        ablation_results.append(result_row)

    # 7. Print and save summary of results
    print("\n--- Ablation Study Summary ---")
    summary_df = pd.DataFrame(ablation_results)
    
    # Define column order for better readability
    main_cols = ['Scenario', 'best_cv_auroc', 'test_auroc', 'test_accuracy', 'test_balanced_accuracy', 'test_f1_score']
    feature_cols = ['Num DNA Features', 'Num RNA Features']
    param_cols = ['best_params']
    
    # Filter out columns that might not exist, then combine
    final_cols = [c for c in main_cols if c in summary_df.columns]
    final_cols += [c for c in feature_cols if c in summary_df.columns]
    final_cols += [c for c in param_cols if c in summary_df.columns]
    
    summary_df = summary_df[final_cols]
    
    print(summary_df.to_string())

    if args.output_csv:
        summary_df.to_csv(args.output_csv, index=False)
        print(f"\nResults saved to {args.output_csv}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Ablation study for NN models with CV and hyperparameter tuning.')
    
    parser.add_argument('--anndata_path', type=str, default='data/combined_data.h5ad', help='Path to the .h5ad data file.')
    parser.add_argument('--test_size', type=float, default=0.2, help='Proportion for the final hold-out test set.')
    parser.add_argument('--output_csv', type=str, default='nn_ablation_results.csv', help='Path to save the output results CSV file.')
    
    parser.add_argument('--model_type', type=str, required=True, choices=['continuous', 'attention'], help='Type of model to train.')
    
    parser.add_argument('--param_grid_json', type=str, required=True, help="Path to a JSON file with the hyperparameter grid.")
    
    parser.add_argument('--k_folds', type=int, default=5, help='Number of folds for cross-validation during tuning.')
    parser.add_argument('--n_epochs', type=int, default=50, help='Number of epochs to train.')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size for training.')

    args = parser.parse_args()
    main(args) 