"""Evaluate saved checkpoints on the test set and report mean ± std across seeds.

Usage:
    python evaluate.py --tag "k10-h64-lr0.001-l3-dp0.3"
    python evaluate.py --tag "k10-h64-lr0.001-l3-dp0.3" --config model_configs/SimpleGNN/default.yaml
"""
import argparse
import os
import yaml
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
from sklearn.metrics import f1_score, roc_auc_score, RocCurveDisplay

from src.utils.config import (PHYLO_DIST_MATRIX, LABELS_RM_METABOLITES,
                               SIMPLE_GNN_XP_RESULTS)
from src.models.SimpleGNN import SimpleGNN
from src.data.build_graph import build_graph

METABOLITE_COLUMNS = ['b12', 'heme', 'folate', 'biotin']

MODEL_REGISTRY = {
    'SimpleGNN': SimpleGNN,
}


def evaluate_test(model, data, mask, metabolite_columns):
    model.eval()
    with torch.no_grad():
        out   = model(data.x, data.edge_index, data.edge_weight, edge_attr=data.edge_attr)
        probs = torch.sigmoid(out[mask]).cpu().numpy()
        preds = (probs > 0.5).astype(int)
        true  = data.y[mask].cpu().numpy().astype(int)

    results = {}
    for i, trait in enumerate(metabolite_columns):
        results[trait] = {
            'f1':    f1_score(true[:, i], preds[:, i], zero_division=0),
            'auroc': roc_auc_score(true[:, i], probs[:, i]),
        }
    results['macro'] = {
        'f1':    np.mean([v['f1']    for v in results.values()]),
        'auroc': np.mean([v['auroc'] for v in results.values()]),
    }
    return results, probs, true


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='model_configs/SimpleGNN/default.yaml')
    parser.add_argument('--tag', required=True,
                        help='Hyperparams portion of run name, e.g. "k10-h64-lr0.001-l3-dp0.3"')
    parser.add_argument('--plot', action='store_true', help='Show ROC curves')
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    device = torch.device('mps' if torch.backends.mps.is_available() else 'cpu')
    print(f"Device: {device}")

    labels      = pd.read_csv(LABELS_RM_METABOLITES)
    dist_matrix = np.load(PHYLO_DIST_MATRIX)

    # Find all seed run folders matching the config tag
    seed_runs = sorted([
        d for d in os.listdir(SIMPLE_GNN_XP_RESULTS)
        if args.tag in d and '-seed' in d
    ])
    if not seed_runs:
        raise FileNotFoundError(
            f"No seed runs found in {SIMPLE_GNN_XP_RESULTS} matching tag '{args.tag}'"
        )
    print(f"\nFound {len(seed_runs)} seed run(s):")
    for r in seed_runs:
        print(f"  {r}")

    all_results = {}
    last_probs, last_true = None, None   # kept for optional ROC plot

    for run_name in seed_runs:
        ckpt = SIMPLE_GNN_XP_RESULTS / run_name / 'model.pt'

        # Rebuild graph with fixed split seed (matches training)
        data = build_graph(
            dist_matrix        = dist_matrix,
            labels_df          = labels,
            metabolite_columns = METABOLITE_COLUMNS,
            k                  = cfg['k'],
            edge_weight        = cfg.get('edge_weight', False),
            edge_attr          = cfg.get('edge_attr', False),
            seed               = 19,
            device             = device,
        )

        model_name = cfg.get('model', 'SimpleGNN')
        ModelClass = MODEL_REGISTRY[model_name]
        model_kwargs = dict(input_dim=cfg['input_dim'], hidden_dim=cfg['hidden_dim'],
                            output_dim=cfg['output_dim'], num_layers=cfg['num_layers'],
                            p_dropout=cfg['p_dropout'])
        if model_name == 'GAT':
            model_kwargs['heads']    = cfg['heads']
            model_kwargs['edge_dim'] = 1 if cfg.get('edge_attr') else None
        model = ModelClass(**model_kwargs).to(device)
        model.load_state_dict(torch.load(ckpt, map_location=device))

        results, probs, true = evaluate_test(model, data, data.test_mask, METABOLITE_COLUMNS)
        all_results[run_name] = results
        last_probs, last_true = probs, true

        pd.DataFrame(results).T.to_csv(SIMPLE_GNN_XP_RESULTS / run_name / 'test_results.csv')

        print(f"\n{run_name}")
        print(f"{'Trait':<10} {'F1':>6} {'AUROC':>7}")
        print("-" * 26)
        for trait, metrics in results.items():
            print(f"{trait:<10} {metrics['f1']:>6.3f} {metrics['auroc']:>7.3f}")

    # Mean ± std — the number you report as the SimpleGNN baseline
    traits = METABOLITE_COLUMNS + ['macro']
    print(f"\n{'=' * 55}")
    print(f"{'Trait':<10} {'AUROC mean':>12} {'± std':>6} {'F1 mean':>9} {'± std':>6}")
    print(f"{'=' * 55}")
    for trait in traits:
        aurocs = [all_results[r][trait]['auroc'] for r in seed_runs]
        f1s    = [all_results[r][trait]['f1']    for r in seed_runs]
        print(f"{trait:<10} {np.mean(aurocs):>12.3f} {np.std(aurocs):>6.3f} "
              f"{np.mean(f1s):>9.3f} {np.std(f1s):>6.3f}")

    if args.plot and last_probs is not None:
        fig, axes = plt.subplots(1, len(METABOLITE_COLUMNS), figsize=(16, 4))
        for i, (trait, ax) in enumerate(zip(METABOLITE_COLUMNS, axes)):
            RocCurveDisplay.from_predictions(last_true[:, i], last_probs[:, i],
                                             ax=ax, name=trait)
            ax.set_title(trait)
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    main()
