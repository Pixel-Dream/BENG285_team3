import argparse
import pandas as pd
from sklearn.linear_model import SGDClassifier
from data_loader import load_and_preprocess_data
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score
import numpy as np

def construct_ablation_feature_matrix(adata, expression_type, mutation_type):
    """
    Constructs a feature matrix for a specific ablation scenario.
    This logic is kept within the training script as requested.
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
        return expr_df.join(mut_df, how='left').fillna(0)
    elif expr_df is not None:
        return expr_df
    elif mut_df is not None:
        return mut_df.reindex(sample_ids).fillna(0)
    else:
        # Return an empty dataframe if no features are selected
        return pd.DataFrame(index=sample_ids)

def main(args):
    """Main training script for Linear SVM ablation study."""
    # 1. Load data
    print("Loading data...")
    adata = load_and_preprocess_data(args.anndata_path)
    labels = adata.obs['survival_3yr_label'].values

    # 2. Define fixed hyperparameters for the model
    # best_params with overfitting
    best_params = {
        'alpha': 0.001, 'l1_ratio': 0.5, 'learning_rate': 'optimal',
        'eta0': 0, 'warm_start': True, 'loss': 'hinge', 'penalty': 'elasticnet',
        'max_iter': 1000, 'tol': 1e-3, 'random_state': 42
    }
    
    # best_params without overfitting
    # best_params = {
    #     'alpha': 0.5, 'l1_ratio': 0.5, 'learning_rate': 'constant',
    #     'eta0': 0.0001, 'warm_start': True, 'loss': 'hinge', 
    #     'penalty': 'elasticnet', 'max_iter': 1000, 'early_stopping': True,
    #     'random_state': 42
    # }

    # 3. Define ablation study scenarios
    ablation_scenarios = [
        # RNA only
        {'name': 'RNA (normalized)', 'expression': 'normalized', 'mutation': None},
        {'name': 'RNA (raw)', 'expression': 'raw', 'mutation': None},
        # Mutations only
        {'name': 'Mutations (full)', 'expression': None, 'mutation': 'jiaming_data_full'},
        {'name': 'Mutations (colinear)', 'expression': None, 'mutation': 'jiaming_data_colinear'},
        {'name': 'Mutations (intersect)', 'expression': None, 'mutation': 'jiaming_data_intersect'},
        # RNA (normalized) + Mutations
        {'name': 'RNA (norm) + Mut (full)', 'expression': 'normalized', 'mutation': 'jiaming_data_full'},
        {'name': 'RNA (norm) + Mut (colinear)', 'expression': 'normalized', 'mutation': 'jiaming_data_colinear'},
        {'name': 'RNA (norm) + Mut (intersect)', 'expression': 'normalized', 'mutation': 'jiaming_data_intersect'},
        # RNA (raw) + Mutations
        {'name': 'RNA (raw) + Mut (full)', 'expression': 'raw', 'mutation': 'jiaming_data_full'},
        {'name': 'RNA (raw) + Mut (colinear)', 'expression': 'raw', 'mutation': 'jiaming_data_colinear'},
        {'name': 'RNA (raw) + Mut (intersect)', 'expression': 'raw', 'mutation': 'jiaming_data_intersect'},
    ]
    
    # 4. Run ablation study
    print("--- Starting SVM Ablation Study ---")
    ablation_results = []

    for scenario in ablation_scenarios:
        print(f"\n--- Running scenario: {scenario['name']} ---")

        feature_matrix = construct_ablation_feature_matrix(
            adata,
            expression_type=scenario['expression'],
            mutation_type=scenario['mutation']
        )
        
        if feature_matrix.empty or feature_matrix.shape[1] == 0:
            print(f"Skipping scenario '{scenario['name']}' due to having 0 features.")
            continue

        print(f"Feature matrix shape: {feature_matrix.shape}")
        
        # Perform a single train-test split
        X_train, X_test, y_train, y_test = train_test_split(
            feature_matrix, labels, test_size=args.test_size, random_state=42, stratify=labels
        )

        # Create a pipeline with scaling and the model
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('model', SGDClassifier(**best_params))
        ])

        # Train the model
        pipeline.fit(X_train, y_train)
        
        # Evaluate on the test set
        y_pred = pipeline.predict(X_test)
        y_score = pipeline.decision_function(X_test)

        # Calculate metrics
        metrics = {
            'accuracy': accuracy_score(y_test, y_pred),
            'balanced_accuracy': balanced_accuracy_score(y_test, y_pred),
            'f1_score': f1_score(y_test, y_pred, average='binary'),
            'auroc': roc_auc_score(y_test, y_score)
        }
        
        # Get feature importance
        importance = np.abs(pipeline.named_steps['model'].coef_[0])
        feature_names = feature_matrix.columns
        sorted_idx = np.argsort(importance)[::-1]
        top_features = feature_names[sorted_idx][:5].tolist()

        metrics['Scenario'] = scenario['name']
        metrics['Num Features'] = feature_matrix.shape[1]
        metrics['Top 5 Features'] = top_features
        ablation_results.append(metrics)

    # 5. Print and save summary of results
    print("\n--- Ablation Study Summary ---")
    summary_df = pd.DataFrame(ablation_results)
    
    column_order = ['Scenario', 'accuracy', 'balanced_accuracy', 'f1_score', 'auroc', 'Num Features', 'Top 5 Features']
    summary_df = summary_df[column_order]
    
    print(summary_df.to_string())

    if args.output_csv:
        summary_df.to_csv(args.output_csv, index=False)
        print(f"\nResults saved to {args.output_csv}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run an ablation study for a Linear SVM model.')
    
    parser.add_argument('--anndata_path', type=str, default='data/combined_data.h5ad', help='Path to the .h5ad data file.')
    parser.add_argument('--test_size', type=float, default=0.2, help='Proportion of the dataset to use for the final hold-out test set.')
    parser.add_argument('--output_csv', type=str, default='svm_ablation_results_overfitted.csv', help='Path to save the output results CSV file.')

    args = parser.parse_args()
    main(args) 