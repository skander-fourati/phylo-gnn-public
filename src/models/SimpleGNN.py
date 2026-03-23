import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv

# num_layers = total layers (first + middles + last), minimum 2
class SimpleGNN(nn.Module):
    def __init__(self,
                 input_dim,
                 hidden_dim,
                 output_dim,
                 num_layers,
                 p_dropout=0.5):

        super().__init__()
        self.dropout = p_dropout  # store so forward can access it

        layers = [GCNConv(input_dim, hidden_dim)]          # first: input → hidden
        for _ in range(num_layers - 2):                    # middles: hidden → hidden (0 if num_layers=2)
            layers.append(GCNConv(hidden_dim, hidden_dim))
        layers.append(GCNConv(hidden_dim, output_dim))     # last: hidden → output
        self.convs = nn.ModuleList(layers)                 # ModuleList registers params with PyTorch

    def forward(self, x, edge_index, edge_weight=None, edge_attr=None):
        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index, edge_weight)
            if i < len(self.convs) - 1:   # relu + dropout between layers, not after last
                x = F.relu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        return x  # raw logits → BCEWithLogitsLoss applies sigmoid