import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv

# num_layers = total layers (first + middles + last), minimum 2
class GAT(nn.Module):
    def __init__(self,
                 input_dim,
                 hidden_dim,
                 output_dim,
                 num_layers,
                 heads,
                 p_dropout=0.5,
                 edge_dim=None):   # set to 1 to feed phylogenetic distances into attention

        super().__init__()
        self.dropout = p_dropout

        # concat=True: stack all head outputs end-to-end → output is [N, hidden_dim * heads]
        layers = [GATConv(input_dim, hidden_dim, heads=heads, concat=True, edge_dim=edge_dim)]
        for _ in range(num_layers - 2):
            layers.append(GATConv(hidden_dim*heads, hidden_dim, heads=heads, concat=True, edge_dim=edge_dim))

        # concat=False: average all head outputs → output is [N, out_channels], used on last layer only
        layers.append(GATConv(hidden_dim*heads, output_dim, heads=heads, concat=False, edge_dim=edge_dim))
        self.convs = nn.ModuleList(layers)

    def forward(self, x, edge_index, edge_weight=None, edge_attr=None):
        # edge_weight ignored (GCNConv concept) — edge_attr feeds into attention if edge_dim was set
        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index, edge_attr=edge_attr)
            if i < len(self.convs) - 1:   # relu + dropout between layers, not after last
                x = F.relu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        return x  # raw logits → BCEWithLogitsLoss applies sigmoid