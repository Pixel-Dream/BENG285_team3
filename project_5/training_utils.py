import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score, f1_score
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterator, *args, **kwargs):
        return iterator
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
import pandas as pd
import joblib
from sklearn.model_selection import ParameterGrid
import itertools


class OmicsDataset(Dataset):
    """Dataset for multi-omics data."""
    def __init__(self, dna_features, rna_features, labels):
        self.dna_features = torch.FloatTensor(dna_features)
        self.rna_features = torch.FloatTensor(rna_features)
        self.labels = torch.LongTensor(labels)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return {
            'dna': self.dna_features[idx],
            'rna': self.rna_features[idx],
            'label': self.labels[idx]
        }

def train_epoch(model, dataloader, criterion, optimizer, device):
    model.train()
    running_loss = 0.0
    all_labels = []
    all_preds = []

    for batch in tqdm(dataloader, desc="Training", leave=False):
        dna = batch['dna'].to(device)
        rna = batch['rna'].to(device)
        labels = batch['label'].to(device)

        optimizer.zero_grad()
        outputs = model(dna, rna)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        running_loss += loss.item() * dna.size(0)
        
        preds = torch.argmax(outputs, dim=1)
        all_labels.extend(labels.cpu().numpy())
        all_preds.extend(preds.cpu().numpy())

    epoch_loss = running_loss / len(dataloader.dataset)
    epoch_acc = accuracy_score(all_labels, all_preds)
    return epoch_loss, epoch_acc

def evaluate_model(model, dataloader, criterion, device):
    model.eval()
    running_loss = 0.0
    all_labels = []
    all_probs = []
    
    with torch.no_grad():
        for batch in tqdm(dataloader, desc="Evaluating", leave=False):
            dna = batch['dna'].to(device)
            rna = batch['rna'].to(device)
            labels = batch['label'].to(device)

            outputs = model(dna, rna)
            loss = criterion(outputs, labels)
            
            running_loss += loss.item() * dna.size(0)
            
            probs = torch.softmax(outputs, dim=1)[:, 1]
            all_labels.extend(labels.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())

    val_loss = running_loss / len(dataloader.dataset)
    val_auc = roc_auc_score(all_labels, all_probs)
    
    return val_loss, val_auc

def k_fold_cross_validation(feature_matrix, labels, model_class, model_kwargs,
                            k=5, n_epochs=50, batch_size=32, learning_rate=0.001):
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)

    fold_results = []
    
    dna_cols = model_kwargs.pop('dna_feature_names')
    rna_cols = model_kwargs.pop('rna_feature_names')

    X_dna = feature_matrix[dna_cols].values
    X_rna = feature_matrix[rna_cols].values
    
    for fold, (train_idx, val_idx) in enumerate(skf.split(feature_matrix, labels)):
        print(f"--- Fold {fold+1}/{k} ---")
        
        X_dna_train, X_dna_val = X_dna[train_idx], X_dna[val_idx]
        X_rna_train, X_rna_val = X_rna[train_idx], X_rna[val_idx]
        y_train, y_val = labels[train_idx], labels[val_idx]

        scaler_dna = StandardScaler().fit(X_dna_train)
        X_dna_train = scaler_dna.transform(X_dna_train)
        X_dna_val = scaler_dna.transform(X_dna_val)

        scaler_rna = StandardScaler().fit(X_rna_train)
        X_rna_train = scaler_rna.transform(X_rna_train)
        X_rna_val = scaler_rna.transform(X_rna_val)
        
        model_kwargs['n_dna_features'] = X_dna_train.shape[1]
        model_kwargs['n_rna_features'] = X_rna_train.shape[1]
        
        model = model_class(**model_kwargs).to(device)
        criterion = nn.CrossEntropyLoss()
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        
        train_dataset = OmicsDataset(X_dna_train, X_rna_train, y_train)
        val_dataset = OmicsDataset(X_dna_val, X_rna_val, y_val)
        
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        
        history = {'train_loss': [], 'train_acc': [], 'val_loss': [], 'val_auc': []}
        
        for epoch in range(n_epochs):
            train_loss, train_acc = train_epoch(model, train_loader, criterion, optimizer, device)
            val_loss, val_auc = evaluate_model(model, val_loader, criterion, device)
            
            history['train_loss'].append(train_loss)
            history['train_acc'].append(train_acc)
            history['val_loss'].append(val_loss)
            history['val_auc'].append(val_auc)
            
            print(f"Epoch {epoch+1}/{n_epochs} - "
                  f"Train Loss: {train_loss:.4f}, Train Acc: {train_acc:.4f}, "
                  f"Val Loss: {val_loss:.4f}, Val AUC: {val_auc:.4f}")

        fold_results.append({
            'history': history,
            'final_auc': val_auc,
            'model': model.state_dict() # Saving model state for potential later use
        })

    all_aucs = [r['final_auc'] for r in fold_results]
    avg_metrics = {
        'avg_auc': np.mean(all_aucs),
        'std_auc': np.std(all_aucs)
    }
    
    return {'fold_results': fold_results, 'average_metrics': avg_metrics}

def plot_cv_results(cv_results):
    """Plots the learning curves for each fold of cross-validation."""
    
    fold_histories = [r['history'] for r in cv_results['fold_results']]
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Cross-Validation Training History', fontsize=16)
    
    # Plot training loss
    for i, history in enumerate(fold_histories):
        sns.lineplot(x=range(len(history['train_loss'])), y=history['train_loss'], 
                     ax=axes[0, 0], label=f"Fold {i+1}")
    axes[0, 0].set_title('Training Loss')
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].set_ylabel('Loss')
    
    # Plot training accuracy
    for i, history in enumerate(fold_histories):
        sns.lineplot(x=range(len(history['train_acc'])), y=history['train_acc'], 
                     ax=axes[0, 1], label=f"Fold {i+1}")
    axes[0, 1].set_title('Training Accuracy')
    axes[0, 1].set_xlabel('Epoch')
    axes[0, 1].set_ylabel('Accuracy')

    # Plot validation loss
    for i, history in enumerate(fold_histories):
        sns.lineplot(x=range(len(history['val_loss'])), y=history['val_loss'], 
                     ax=axes[1, 0], label=f"Fold {i+1}")
    axes[1, 0].set_title('Validation Loss')
    axes[1, 0].set_xlabel('Epoch')
    axes[1, 0].set_ylabel('Loss')

    # Plot validation AUC
    for i, history in enumerate(fold_histories):
        sns.lineplot(x=range(len(history['val_auc'])), y=history['val_auc'], 
                     ax=axes[1, 1], label=f"Fold {i+1}")
    axes[1, 1].set_title('Validation AUC')
    axes[1, 1].set_xlabel('Epoch')
    axes[1, 1].set_ylabel('AUC')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

def run_sklearn_grid_search(
    X: pd.DataFrame, 
    y: np.ndarray, 
    model, 
    param_grid: dict, 
    test_size: float = 0.2, 
    cv: int = 5,
    random_state: int = 42
    ):
    """
    Runs GridSearchCV for a scikit-learn model, including a final evaluation on a hold-out test set.
    """
    # Split data into training and hold-out test sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )

    print(f"Training set size: {X_train.shape[0]}")
    print(f"Test set size: {X_test.shape[0]}")

    # Create a pipeline with scaling and the model
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('model', model)
    ])
    
    # Adjust param_grid to match pipeline format (e.g., 'model__alpha')
    pipeline_param_grid = {f'model__{key}': value for key, value in param_grid.items()}

    # Perform grid search with cross-validation on the training set
    print(f"Running GridSearchCV with {cv}-fold cross-validation...")
    grid_search = GridSearchCV(
        pipeline,
        param_grid=pipeline_param_grid,
        cv=cv,
        scoring='accuracy',
        n_jobs=-1,
        verbose=1
    )
    grid_search.fit(X_train, y_train)

    # Evaluate the best model on the hold-out test set
    best_model = grid_search.best_estimator_
    y_pred = best_model.predict(X_test)
    
    # For AUROC, we need a score, not a prediction. Use decision_function for SVMs.
    if hasattr(best_model, "decision_function"):
        y_score = best_model.decision_function(X_test)
    else:
        # Fallback for models without decision_function (though SGDClassifier has it)
        y_score = y_pred

    test_metrics = {
        'accuracy': accuracy_score(y_test, y_pred),
        'balanced_accuracy': balanced_accuracy_score(y_test, y_pred),
        'f1_score': f1_score(y_test, y_pred, average='binary'),
        'auroc': roc_auc_score(y_test, y_score)
    }

    # Extract feature importance if available
    feature_importance = None
    if hasattr(best_model.named_steps['model'], 'coef_'):
        # Coefficients are the importance for linear models
        importance = np.abs(best_model.named_steps['model'].coef_[0])
        feature_names = X.columns
        
        sorted_idx = np.argsort(importance)[::-1]
        
        feature_importance = pd.DataFrame({
            'feature': feature_names[sorted_idx],
            'importance': importance[sorted_idx]
        })

    results = {
        'best_params': {k.replace('model__', ''): v for k, v in grid_search.best_params_.items()},
        'best_cv_score': grid_search.best_score_,
        'test_metrics': test_metrics,
        'best_model': best_model,
        'feature_importance': feature_importance
    }

    print("\n--- Grid Search Complete ---")
    print(f"Best parameters found: {results['best_params']}")
    print(f"Best cross-validation accuracy: {results['best_cv_score']:.4f}")
    print("\nHold-out Test Set Performance:")
    for metric, value in results['test_metrics'].items():
        print(f"- {metric.replace('_', ' ').title()}: {value:.4f}")
    
    if feature_importance is not None:
        print("\nTop 5 most important features:")
        print(feature_importance.head(5))

    return results 

def evaluate_sklearn_model(
    X: pd.DataFrame,
    y: np.ndarray,
    model,
    test_size: float = 0.2,
    random_state: int = 42
    ):
    """
    Trains and evaluates a scikit-learn model on a hold-out test set without cross-validation.
    """
    # Split data into training and hold-out test sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )

    print(f"Training set size: {X_train.shape[0]}")
    print(f"Test set size: {X_test.shape[0]}")

    # Create a pipeline with scaling and the model
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('model', model)
    ])

    # Train the model on the full training set
    print("Training model...")
    pipeline.fit(X_train, y_train)

    # Evaluate the best model on the hold-out test set
    y_pred = pipeline.predict(X_test)

    if hasattr(pipeline, "decision_function"):
        y_score = pipeline.decision_function(X_test)
    else:
        y_score = y_pred

    test_metrics = {
        'accuracy': accuracy_score(y_test, y_pred),
        'balanced_accuracy': balanced_accuracy_score(y_test, y_pred),
        'f1_score': f1_score(y_test, y_pred, average='binary'),
        'auroc': roc_auc_score(y_test, y_score)
    }

    # Extract feature importance if available
    feature_importance = None
    if hasattr(pipeline.named_steps['model'], 'coef_'):
        importance = np.abs(pipeline.named_steps['model'].coef_[0])
        feature_names = X.columns
        
        sorted_idx = np.argsort(importance)[::-1]
        
        feature_importance = pd.DataFrame({
            'feature': feature_names[sorted_idx],
            'importance': importance[sorted_idx]
        })

    results = {
        'test_metrics': test_metrics,
        'feature_importance': feature_importance
    }

    print("\nHold-out Test Set Performance:")
    for metric, value in results['test_metrics'].items():
        print(f"- {metric.replace('_', ' ').title()}: {value:.4f}")

    if feature_importance is not None:
        print("\nTop 5 most important features:")
        print(feature_importance.head(5))

    return results

def tune_pytorch_hyperparameters(
    feature_matrix: pd.DataFrame, 
    labels: np.ndarray,
    model_class,
    param_grid: dict,
    static_model_kwargs: dict,
    k: int = 5,
    n_epochs: int = 50,
    batch_size: int = 32
    ):
    """
    Performs grid search hyperparameter tuning for a PyTorch model using k-fold cross-validation.
    """
    grid = ParameterGrid(param_grid)
    best_score = -1
    best_params = None
    all_results = []

    print(f"Starting hyperparameter tuning with {len(list(grid))} combinations...")

    for params in grid:
        print(f"\n--- Testing parameters: {params} ---")
        
        # Combine static kwargs with current grid params
        current_model_kwargs = static_model_kwargs.copy()
        
        # Separate learning rate as it's not a model kwarg
        learning_rate = params.pop('learning_rate', 0.001)
        
        # The rest of params are for the model
        current_model_kwargs.update(params)

        # Run k-fold CV with the current parameter set
        cv_results = k_fold_cross_validation(
            feature_matrix=feature_matrix.copy(),
            labels=labels,
            model_class=model_class,
            model_kwargs=current_model_kwargs,
            k=k,
            n_epochs=n_epochs,
            batch_size=batch_size,
            learning_rate=learning_rate
        )

        avg_auc = cv_results['average_metrics']['avg_auc']
        
        # Add back learning rate for logging
        params['learning_rate'] = learning_rate
        
        all_results.append({
            'params': params,
            'avg_auc': avg_auc,
            'std_auc': cv_results['average_metrics']['std_auc']
        })

        print(f"Parameters: {params} -> Avg. CV AUC: {avg_auc:.4f}")

        if avg_auc > best_score:
            best_score = avg_auc
            best_params = params

    print("\n--- Hyperparameter Tuning Complete ---")
    print(f"Best parameters found: {best_params}")
    print(f"Best average cross-validation AUC: {best_score:.4f}")

    return {'best_params': best_params, 'best_score': best_score, 'all_results': all_results} 