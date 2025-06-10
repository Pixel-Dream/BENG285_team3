import argparse
import pandas as pd
from model.models import ContinuousFeatureNN, AttentionFeatureNN
from data_loader import load_and_preprocess_data, construct_feature_matrix, get_feature_counts
from training_utils import k_fold_cross_validation, plot_cv_results, tune_pytorch_hyperparameters
import json

def main(args):
    """Main training script."""
    # 1. Load and preprocess data
    print("Loading data...")
    adata = load_and_preprocess_data(args.data_path)
    
    # 2. Construct feature matrix
    print("Constructing feature matrix...")
    feature_matrix = construct_feature_matrix(
        adata,
        expression_type=args.expression_type,
        mutation_data_source=args.mutation_data_source,
        mutation_data_identifier=args.mutation_data_identifier
    )
    
    labels = adata.obs['survival_3yr_label'].values
    
    # 3. Get feature counts and define static model args
    feature_counts = get_feature_counts(feature_matrix)
    print(f"Total features: {feature_counts['total_features']}")
    print(f"RNA features: {feature_counts['n_rna_features']}")
    print(f"DNA features: {feature_counts['n_dna_features']}")

    static_model_kwargs = {
        'dna_feature_names': feature_counts['dna_feature_names'],
        'rna_feature_names': feature_counts['rna_feature_names']
    }

    # 4. Model selection
    if args.model_type == 'continuous':
        model_class = ContinuousFeatureNN
    elif args.model_type == 'attention':
        model_class = AttentionFeatureNN
    else:
        raise ValueError(f"Unknown model type: {args.model_type}")

    # 5. Run tuning or a single training run
    if args.tune:
        if not args.param_grid_json:
            raise ValueError("Must provide --param_grid_json when --tune is set.")
        with open(args.param_grid_json, 'r') as f:
            param_grid = json.load(f)
        print("Starting hyperparameter tuning...")
        
        tune_results = tune_pytorch_hyperparameters(
            feature_matrix=feature_matrix,
            labels=labels,
            model_class=model_class,
            param_grid=param_grid,
            static_model_kwargs=static_model_kwargs,
            k=args.k_folds,
            n_epochs=args.n_epochs,
            batch_size=args.batch_size
        )
        print("\nTuning summary:")
        print(pd.DataFrame(tune_results['all_results']).sort_values('avg_auc', ascending=False))

    else:
        print(f"Starting a single {args.k_folds}-fold cross-validation run...")
        model_kwargs = {
            **static_model_kwargs,
            'hidden_dim': args.hidden_dim,
            'dropout_rate': args.dropout_rate,
        }
        
        cv_results = k_fold_cross_validation(
            feature_matrix=feature_matrix,
            labels=labels,
            model_class=model_class,
            model_kwargs=model_kwargs,
            k=args.k_folds,
            n_epochs=args.n_epochs,
            batch_size=args.batch_size,
            learning_rate=args.learning_rate
        )
        
        avg_metrics = cv_results['average_metrics']
        print(f"\nAverage {args.k_folds}-Fold CV Performance:")
        print(f"AUC: {avg_metrics['avg_auc']:.3f} ± {avg_metrics['std_auc']:.3f}")
        
        if args.plot_results:
            print("Plotting results...")
            plot_cv_results(cv_results)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train a multi-omics neural network.')
    
    # Data arguments
    parser.add_argument('--data_path', type=str, default='data/combined_data.h5ad',
                        help='Path to the .h5ad data file.')
    parser.add_argument('--expression_type', type=str, default='normalized',
                        choices=['normalized', 'raw'],
                        help='Type of expression data to use.')
    parser.add_argument('--mutation_data_source', type=str, default='anndata',
                        choices=['anndata', 'csv', 'none'],
                        help="Source for mutation data ('anndata', 'csv', or 'none').")
    parser.add_argument('--mutation_data_identifier', type=str, default='jiaming_data_full',
                        help='Identifier for mutation data. For anndata, key in .uns. For csv, path to file.')

    # Model arguments
    parser.add_argument('--model_type', type=str, default='attention',
                        choices=['continuous', 'attention'],
                        help='Type of model to train.')
    # These are now for single runs, tuning uses JSON
    parser.add_argument('--hidden_dim', type=int, default=128,
                        help='Hidden dimension size in the model (for single run).')
    parser.add_argument('--dropout_rate', type=float, default=0.3,
                        help='Dropout rate for the model (for single run).')
    parser.add_argument('--learning_rate', type=float, default=0.001,
                        help='Learning rate for the optimizer (for single run).')

    # Training arguments
    parser.add_argument('--k_folds', type=int, default=5,
                        help='Number of folds for cross-validation.')
    parser.add_argument('--n_epochs', type=int, default=50,
                        help='Number of epochs to train.')
    parser.add_argument('--batch_size', type=int, default=32,
                        help='Batch size for training.')

    # Tuning arguments
    parser.add_argument('--tune', action='store_true',
                        help='If set, run hyperparameter tuning instead of a single training run.')
    parser.add_argument('--param_grid_json', type=str, default=None,
                        help="Path to a JSON file containing the hyperparameter grid for tuning.")

    # Other arguments
    parser.add_argument('--plot_results', action='store_true',
                        help='If set, plot the cross-validation results (for single run).')

    args = parser.parse_args()
    
    if args.mutation_data_source.lower() == 'none':
        args.mutation_data_source = None
        
    main(args)
