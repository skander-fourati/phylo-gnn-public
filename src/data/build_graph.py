import numpy as np
import torch
from torch_geometric.data import Data
from sklearn.model_selection import train_test_split


def build_graph(dist_matrix, labels_df, metabolite_columns, k,
                edge_weight=False, edge_attr=False, seed=19, device='cpu'):
    """Build a PyG Data object from a phylogenetic distance matrix.

    Constructs a k-NN graph where each node is a species, edges connect the
    k nearest phylogenetic neighbours, and node features are binary metabolic
    trait labels. Val/test node features are zeroed out (semi-supervised setup).

    Args:
        dist_matrix:        (N x N) numpy array of pairwise phylogenetic distances.
        labels_df:          DataFrame with one row per species; must contain metabolite_columns.
        metabolite_columns: List of column names used as node features and labels.
        k:                  Number of nearest neighbours per node.
        edge_weight:        If True, edge weights = 1 / phylogenetic_distance.
        seed:               Random state for train/val/test split (default 19).
        device:             torch device string or object.

    Returns:
        PyG Data object with x, y, edge_index, edge_weight (or None),
        train_mask, val_mask, test_mask — all moved to device.
    """
    # Build k-NN edges from distance matrix
    edges, distances = [], []
    for i in range(len(dist_matrix)):
        top_k = np.argsort(dist_matrix[i])[1:k + 1]   # skip self (index 0)
        for j in top_k:
            edges.append([i, j])
            distances.append(dist_matrix[i][j])

    edge_index = torch.tensor(np.array(edges).T, dtype=torch.long)

    # Node features = labels (val/test nodes get zeroed out below)
    x = torch.tensor(labels_df[metabolite_columns].values, dtype=torch.float)
    y = x.clone()

    data = Data(x=x, edge_index=edge_index, y=y)

    # Train 60 / val 20 / test 20 split
    indices = np.arange(len(labels_df))
    train_idx, temp   = train_test_split(indices, test_size=0.4, random_state=seed)
    val_idx, test_idx = train_test_split(temp,    test_size=0.5, random_state=seed)

    train_mask = torch.zeros(data.num_nodes, dtype=torch.bool)
    val_mask   = torch.zeros(data.num_nodes, dtype=torch.bool)
    test_mask  = torch.zeros(data.num_nodes, dtype=torch.bool)
    train_mask[train_idx] = True
    val_mask[val_idx]     = True
    test_mask[test_idx]   = True

    data.train_mask = train_mask
    data.val_mask   = val_mask
    data.test_mask  = test_mask

    # Semi-supervised: hide val and test features so the model can't cheat
    data.x[val_mask | test_mask] = 0

    inv_distances = [1 / d for d in distances]

    # GCNConv: scalar edge weight directly scales aggregated messages
    data.edge_weight = torch.tensor(inv_distances, dtype=torch.float) if edge_weight else None

    # GATConv: edge features fed into attention coefficient computation [num_edges, 1]
    data.edge_attr = torch.tensor(inv_distances, dtype=torch.float).unsqueeze(-1) if edge_attr else None

    # Move everything to device
    data.x          = data.x.to(device)
    data.y          = data.y.to(device)
    data.edge_index = data.edge_index.to(device)
    data.train_mask = data.train_mask.to(device)
    data.val_mask   = data.val_mask.to(device)
    data.test_mask  = data.test_mask.to(device)
    if data.edge_weight is not None:
        data.edge_weight = data.edge_weight.to(device)
    if data.edge_attr is not None:
        data.edge_attr = data.edge_attr.to(device)

    return data
