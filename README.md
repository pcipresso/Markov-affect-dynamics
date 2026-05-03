# Markov-affect-dynamics

Reproducible R workflow for analyzing affect dynamics from intensive longitudinal data using observed-state discrete-time Markov chains.

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
