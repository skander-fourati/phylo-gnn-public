# Phylogeny-Aware GNN for Bacterial Metabolism Prediction

> **Work in progress** — active development as of March 2026.

Predicting bacterial metabolic capabilities from genomic data using graph neural networks that incorporate phylogenetic relationships. Applied to 4,744 human gut bacterial species from the [UHGG v2.0](https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0) dataset.

---

## Motivation

Bacteria in the human gut produce vitamins and cofactors (B12, folate, biotin, heme) that humans cannot synthesize independently. Knowing which species produce which metabolites is critical for understanding microbiome contributions to host health — but experimentally characterizing thousands of species is infeasible.

This project asks: **can we predict these traits computationally, and does knowing the evolutionary history of a species improve those predictions?**

A second question motivates the model design: different traits have very different evolutionary dynamics. Heme biosynthesis is nearly universal and vertically inherited, while restriction-modification (R-M) systems are frequently horizontally transferred between species. A phylogeny-aware model should perform better on conserved traits and worse on HGT-prone ones — testing this is a core part of the analysis.

---

## Dataset

- **Species**: 4,744 UHGG v2.0 species-level representative genomes
- **Phylogeny**: GTDB bac120 IQ-TREE phylogenetic tree → pairwise distance matrix (4,744 × 4,744)
- **Labels** (binary, per species):
  - B12 (cobalamin) biosynthesis
  - Heme biosynthesis
  - Folate biosynthesis
  - Biotin biosynthesis
  - Restriction-modification (R-M) systems

Labels are derived from pre-computed UHGG KEGG module completeness scores (≥75% module completeness = producer). Multi-pathway traits (e.g., B12 has 3 alternative pathways) use max completeness across pathways. Metabolic labels are cross-validated against AGORA2.

---

## Approach

The phylogenetic distance matrix is converted to a k-NN graph (k=10 nearest phylogenetic neighbors). Each node represents a bacterial species; edges connect phylogenetically similar species weighted by inverse distance. GNNs perform semi-supervised node classification: training labels are propagated through the graph to predict held-out species.

**Models compared:**
| Model | Description |
|-------|-------------|
| Logistic Regression / MLP | Non-graph baseline (same features, no phylogeny) |
| SimpleGNN | 3-layer GCN with phylogenetic edge weights |
| GAT | Graph Attention Networks with multi-head attention |
| Pretrained genomic model | Embeddings from a pretrained genomic foundation model |

---

## Preliminary Results

SimpleGNN baseline (macro-averaged across 4 metabolic traits, 3-seed average):

| Trait | AUROC | F1 |
|-------|-------|----|
| Heme | 0.963 | 0.852 |
| Biotin | 0.959 | 0.813 |
| Folate | 0.926 | 0.824 |
| B12 | 0.905 | 0.765 |
| **Macro** | **0.938** | **0.814** |

The trait-level pattern is consistent with evolutionary priors: heme (most conserved) is most predictable, B12 (moderate HGT) is hardest. GAT evaluation and full model comparisons are in progress.

---

## Tech Stack

- **PyTorch Geometric** — graph neural networks
- **Biopython** — phylogenetic tree parsing and distance computation
- **KEGG REST API** — module definitions and completeness parsing
- **GTDB / UHGG** — reference phylogeny and genome catalog

---

## Project Structure

```
src/
  models/        # SimpleGNN, GAT implementations
  data/          # UHGG parser, graph construction, KEGG client
  utils/         # config, HTTP utilities
notebooks/
  pipeline/      # Data processing (labels, phylo tree, R-M systems)
  experiments/   # Training and evaluation notebooks
train.py         # Grid search training runner
evaluate.py      # Test set evaluation with seed averaging
data/processed/  # labels.csv, phylo_distance_matrix.npy
```

---

## Status

| Component | Status |
|-----------|--------|
| Data pipeline | Complete |
| SimpleGNN baseline | Complete |
| GAT model | Trained, evaluation pending |
| Non-graph baseline | In progress |
| Pretrained model comparison | Planned |
| Phylogeny ablation study | Planned |
| Full write-up | Planned |
