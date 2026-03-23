"""Train SimpleGNN (or any registered model) on the phylogenetic GNN task.

Usage:
    # Single run with config defaults
    python train.py --config model_configs/SimpleGNN/default.yaml

    # Single run with a fixed seed (for seed sweep after grid search)
    python train.py --config model_configs/SimpleGNN/default.yaml --seed 0

    # Grid search (overrides grid_search flag in YAML)
    python train.py --config model_configs/SimpleGNN/default.yaml --grid-search
"""
import argparse
import datetime
import yaml
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from itertools import product
from sklearn.metrics import f1_score, roc_auc_score
from torch.utils.tensorboard import SummaryWriter

from src.utils.config import (PHYLO_DIST_MATRIX, LABELS_RM_METABOLITES,
                               RUNS_DIR, MODEL_CONFIGS_DIR)
from src.models.SimpleGNN import SimpleGNN
from src.models.GAT import GAT
from src.data.build_graph import build_graph

METABOLITE_COLUMNS = ['b12', 'heme', 'folate', 'biotin']

# Registry: add new architectures here as the project grows
MODEL_REGISTRY = {
    'SimpleGNN': SimpleGNN,
    'GAT':       GAT,
}

PARAM_GRID = {
    'hidden_dim': [32, 64, 128],
    'num_layers': [2, 3],
    'lr':         [0.01, 0.001],
    'p_dropout':  [0.3, 0.5],
    'heads':      [2,3,4]
}


def evaluate_epoch(model, data, mask, epoch, labels_str, metadata_headers, criterion, writer):
    model.eval()
    with torch.no_grad():
        out      = model(data.x, data.edge_index, data.edge_weight, edge_attr=data.edge_attr)
        val_loss = criterion(out[mask], data.y[mask]).item()
        probs    = torch.sigmoid(out[mask]).cpu().numpy()

    preds = (probs > 0.5).astype(int)
    true  = data.y[mask].cpu().numpy().astype(int)

    results = {'val_loss': val_loss}
    writer.add_scalar('Loss/val', val_loss, epoch)

    for i, trait in enumerate(METABOLITE_COLUMNS):
        auc = roc_auc_score(true[:, i], probs[:, i])
        f1  = f1_score(true[:, i], preds[:, i], zero_division=0)
        writer.add_scalar(f'F1/{trait}',    f1,  epoch)
        writer.add_scalar(f'AUROC/{trait}', auc, epoch)
        results[f'F1/{trait}']    = f1
        results[f'AUROC/{trait}'] = auc

    if epoch % 25 == 0:
        if data.edge_attr is not None:
            embedding = model.convs[0](data.x, data.edge_index, edge_attr=data.edge_attr).cpu()
        else:
            embedding = model.convs[0](data.x, data.edge_index, data.edge_weight).cpu()
        writer.add_embedding(embedding, metadata=labels_str,
                             metadata_header=metadata_headers,
                             global_step=epoch, tag=f'epoch_{epoch}')
    return results


def train_run(cfg, data, labels, device, today,
              hidden_dim=None, num_layers=None, lr=None, p_dropout=None,
              heads=None, seed=None):
    hidden_dim = hidden_dim or cfg['hidden_dim']
    num_layers = num_layers or cfg['num_layers']
    lr         = lr         or cfg['lr']
    p_dropout  = p_dropout  or cfg['p_dropout']
    k          = cfg['k']
    epochs     = cfg['epochs']

    if seed is not None:
        torch.manual_seed(seed)

    model_name = cfg.get('model', 'SimpleGNN')
    ModelClass = MODEL_REGISTRY[model_name]
    model_kwargs = dict(input_dim=cfg['input_dim'], hidden_dim=hidden_dim,
                        output_dim=cfg['output_dim'], num_layers=num_layers,
                        p_dropout=p_dropout)
    if model_name == 'GAT':
        heads = heads or cfg['heads']
        model_kwargs['heads']    = heads
        model_kwargs['edge_dim'] = 1 if cfg.get('edge_attr') else None
    model = ModelClass(**model_kwargs).to(device)

    pos_counts = torch.tensor(labels[METABOLITE_COLUMNS].sum().values, dtype=torch.float)
    neg_counts = len(labels) - pos_counts
    pos_weight = (neg_counts / pos_counts).to(device) if cfg.get('pos_weight') else None

    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    train_losses, val_losses = [], []
    auroc_lists = {trait: [] for trait in METABOLITE_COLUMNS}
    f1_lists    = {trait: [] for trait in METABOLITE_COLUMNS}

    seed_tag  = f'-seed{seed}' if seed is not None else ''
    heads_tag = f'-heads{heads}' if model_name == 'GAT' else ''
    run_name  = f"{today.strftime('%Y-%m-%d_%H%M%S')}-k{k}-h{hidden_dim}-lr{lr}-l{num_layers}-dp{p_dropout}{heads_tag}{seed_tag}"
    results_dir = RUNS_DIR / model_name

    metadata_headers = METABOLITE_COLUMNS + ['genome_id']
    labels_str = labels[metadata_headers].astype(str).values.tolist()

    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / run_name).mkdir(parents=True, exist_ok=True)

    # Save the actual hyperparams used (grid search overrides cfg defaults)
    run_cfg = {**cfg, 'hidden_dim': hidden_dim, 'num_layers': num_layers,
               'lr': lr, 'p_dropout': p_dropout}
    if model_name == 'GAT':
        run_cfg['heads'] = heads
    with open(results_dir / run_name / 'config.yaml', 'w') as f:
        yaml.dump(run_cfg, f)

    writer = SummaryWriter(results_dir / run_name)
    writer.add_graph(model=model, input_to_model=(data.x, data.edge_index))

    for epoch in range(epochs):
        model.train()
        out = model(data.x, data.edge_index, data.edge_weight, edge_attr=data.edge_attr)
        optimizer.zero_grad()
        loss = criterion(out[data.train_mask], data.y[data.train_mask])
        loss.backward()
        optimizer.step()
        writer.add_scalar('Loss/train', loss.item(), epoch)
        train_losses.append(loss.item())

        eval_results = evaluate_epoch(model, data, data.val_mask, epoch,
                                      labels_str, metadata_headers, criterion, writer)
        val_losses.append(eval_results['val_loss'])
        for trait in METABOLITE_COLUMNS:
            auroc_lists[trait].append(eval_results[f'AUROC/{trait}'])
            f1_lists[trait].append(eval_results[f'F1/{trait}'])

        if epoch % 10 == 0:
            print(f"Epoch {epoch:03d} | Train: {loss.item():.4f} | Val: {eval_results['val_loss']:.4f}")

    writer.close()

    save_path = results_dir / run_name / 'model.pt'
    torch.save(model.state_dict(), save_path)
    print(f"Saved → {save_path}")

    best_aurocs = {trait: max(auroc_lists[trait]) for trait in METABOLITE_COLUMNS}
    best_f1s    = {trait: max(f1_lists[trait])    for trait in METABOLITE_COLUMNS}

    hparam_writer = SummaryWriter(results_dir / run_name / 'hparams')
    hparams = {'k': k, 'hidden_dim': hidden_dim, 'lr': lr, 'num_layers': num_layers,
               'pos_weight': cfg.get('pos_weight', False)}
    if model_name == 'GAT':
        hparams['heads']     = heads
        hparams['edge_attr'] = cfg.get('edge_attr', False)
    else:
        hparams['edge_weight'] = cfg.get('edge_weight', False)

    hparam_writer.add_hparams(
        hparams,
        {'best_train_loss': min(train_losses),
         'best_val_loss':   min(val_losses),
         'auroc_b12':       best_aurocs['b12'],
         'auroc_heme':      best_aurocs['heme'],
         'auroc_folate':    best_aurocs['folate'],
         'auroc_biotin':    best_aurocs['biotin'],
         'auroc_macro':     sum(best_aurocs.values()) / len(METABOLITE_COLUMNS)}
    )
    hparam_writer.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='model_configs/SimpleGNN/default.yaml')
    parser.add_argument('--grid-search', action='store_true')
    parser.add_argument('--seed', type=int, default=None,
                        help='Fixed seed for weight init + data split. Use for seed sweep after grid search.')
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    device = torch.device('mps' if torch.backends.mps.is_available() else 'cpu')
    print(f"Device: {device}")

    labels      = pd.read_csv(LABELS_RM_METABOLITES)
    dist_matrix = np.load(PHYLO_DIST_MATRIX)
    today       = datetime.datetime.today()

    data = build_graph(
        dist_matrix        = dist_matrix,
        labels_df          = labels,
        metabolite_columns = METABOLITE_COLUMNS,
        k                  = cfg['k'],
        edge_weight        = cfg.get('edge_weight', False),
        edge_attr          = cfg.get('edge_attr', False),
        seed               = args.seed or 19,
        device             = device,
    )

    run_grid_search = args.grid_search or cfg.get('grid_search', False)

    if run_grid_search:
        keys = list(PARAM_GRID.keys())
        for values in product(*PARAM_GRID.values()):
            params = dict(zip(keys, values))
            train_run(cfg, data, labels, device, today, **params)
    else:
        train_run(cfg, data, labels, device, today, seed=args.seed)


if __name__ == '__main__':
    main()
