# Markov-affect-dynamics

Reproducible R workflow for analyzing affect dynamics from intensive longitudinal data using observed-state discrete-time Markov chains.

**Authors:** Pietro Cipresso & Francesca Borghesi

## Overview

This repository provides a complete and modular pipeline for modeling affect dynamics through discrete-time Markov chains based on observed states.

The workflow connects:

- Raw intensive longitudinal data (e.g., ESM, PANAS)
- State construction and discretization
- Transition count matrix (C matrix)
- Transition probability matrix (P matrix)
- Extraction of Markov-based features
- Visualization (heatmaps, transition graphs)

The implementation follows the methodological framework described in the accompanying tutorial.

## Repository structure

- `Preprocessing_row_data_to_C_matrix.R`  
  Transforms raw time-indexed data into:
  - discrete state sequence
  - transition count matrix (C)

- `Markov_affect_dynamics.R`  
  Takes the transition matrix as input and computes:
  - transition probabilities (P)
  - persistence, entropy, asymmetry, attractor strength
  - graphical outputs (heatmaps, network plots)

## Input data requirements

The workflow assumes a long-format dataset with at least:

- `id` → subject identifier  
- `time` → ordered time index  
- affect measures (e.g., PANAS items, valence/arousal, or derived scores)

Data can come from:

- experience sampling (ESM)
- laboratory paradigms
- psychophysiological recordings
- virtual reality environments

## How to run the workflow

### Step 1 — Preprocessing

Run:

```r
source("Preprocessing_row_data_to_C_matrix.R")
```

This step:

- maps raw affect data into discrete states
- constructs the transition count matrix (C)

---

### Step 2 — Markov analysis

Run:

```r
source("Markov_affect_dynamics.R")
```

This step:

- estimates transition probabilities (P)
- computes Markov-based indices
- generates visualizations

---

## Output

The workflow produces:

- Transition count matrix (C)
- Transition probability matrix (P)
- State-level indices:
  - persistence
  - mobility
  - entropy
  - directional asymmetry
  - attractor strength
- System-level summaries
- Visual outputs:
  - heatmaps
  - transition graphs

---

## Notes

- The workflow is modular: users can enter at the preprocessing stage or directly from an existing transition matrix.
- State construction must be explicitly defined and justified.
- The same pipeline applies across different data modalities once a discrete state mapping is specified.

---

## Reproducibility

The full pipeline can be executed end-to-end using the provided scripts, allowing researchers to move from raw intensive longitudinal data to interpretable transition-based features with minimal manual intervention.

---

## Related publication

Full methodological details are provided in the associated tutorial:

*Discrete-Time Markov Chains for Affect Dynamics: A Reproducible Workflow* (under review)
