import torch
import torch.nn as nn

class ContinuousFeatureNN(nn.Module):
    """Neural network for continuous gene expression and mutation count features."""
    
    def __init__(self, n_dna_features, n_rna_features, hidden_dim=128, dropout_rate=0.3):
        super().__init__()
        
        self.dna_encoder = nn.Sequential(
            nn.Linear(n_dna_features, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU()
        )
        
        self.rna_encoder = nn.Sequential(
            nn.Linear(n_rna_features, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU()
        )
        
        fusion_dim = hidden_dim
        self.fusion = nn.Sequential(
            nn.Linear(fusion_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout_rate)
        )
        
        self.classifier = nn.Linear(hidden_dim // 2, 2)
        
        self._init_weights()
    
    def _init_weights(self):
        """Initialize weights for linear layers."""
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight)
                if module.bias is not None:
                    nn.init.constant_(module.bias, 0)
            elif isinstance(module, nn.BatchNorm1d):
                nn.init.constant_(module.weight, 1)
                nn.init.constant_(module.bias, 0)
    
    def forward(self, dna_features, rna_features):
        """
        Args:
            dna_features: Tensor of shape (batch_size, n_dna_features)
            rna_features: Tensor of shape (batch_size, n_rna_features)
        """
        dna_encoded = self.dna_encoder(dna_features)
        rna_encoded = self.rna_encoder(rna_features)
        
        combined = torch.cat([dna_encoded, rna_encoded], dim=1)
        
        fused = self.fusion(combined)
        output = self.classifier(fused)
        
        return output


class AttentionFeatureNN(nn.Module):
    """Neural network with attention mechanism for feature importance."""
    
    def __init__(self, n_dna_features, n_rna_features, hidden_dim=128, dropout_rate=0.3):
        super().__init__()
        
        self.n_dna_features = n_dna_features
        self.n_rna_features = n_rna_features
        
        self.dna_attention = nn.Sequential(
            nn.Linear(n_dna_features, n_dna_features),
            nn.Sigmoid()
        )
        
        self.rna_attention = nn.Sequential(
            nn.Linear(n_rna_features, n_rna_features),
            nn.Sigmoid()
        )
        
        self.main_network = ContinuousFeatureNN(
            n_dna_features, n_rna_features, hidden_dim, dropout_rate
        )
    
    def forward(self, dna_features, rna_features):
        """Apply attention then process with main network."""
        dna_weights = self.dna_attention(dna_features)
        rna_weights = self.rna_attention(rna_features)
        
        dna_attended = dna_features * dna_weights
        rna_attended = rna_features * rna_weights
        
        return self.main_network(dna_attended, rna_attended)
    
    def get_feature_importance(self, dna_features, rna_features):
        """Get attention weights for feature importance analysis."""
        with torch.no_grad():
            dna_weights = self.dna_attention(dna_features)
            rna_weights = self.rna_attention(rna_features)
            
            return {
                'dna_importance': dna_weights.mean(dim=0).cpu().numpy(),
                'rna_importance': rna_weights.mean(dim=0).cpu().numpy()
            } 