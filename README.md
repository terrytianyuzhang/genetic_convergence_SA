# XConTet for Genetic Convergence Analysis 

## Overview

We developed XConTest, a cross-fitted convergence testing framework designed to evaluate the similarity of transcriptional impacts induced by genetic perturbations. The method is demonstrated using multiple Perturb-seq datasets.

## Repository Structure

Each folder corresponds to a specific section of the paper:

 - `yao_2023/` Application of our method to the dataset reported in Yao et al., Scalable genetic screening for regulatory circuits using compressed Perturb-seq.
- `lalli_2020/` Application to Lalli et al., High-throughput single-cell functional elucidation of neurodevelopmental disease–associated genes reveals convergent mechanisms altering neuronal differentiation.
- `simulation/` Simulation studies evaluating type-I error control and statistical power.
- `converge_by_module/` Simulation studies demonstrating gene-module–level convergence using the XConTest procedure.

Within each folder, a Quarto document `scripts/reproduce_results.qmd` reproduces all corresponding analysis and figures.

## Software Package

The software implementation accompanying this paper is available [here](https://github.com/terrytianyuzhang/HMC)

A vignette illustrating the use of XConTest can be found  [here](https://terrytianyuzhang.github.io/HMC/HMC_convergence.html) 

## Citation

If you use this method or repository in your work, please cite:

**Zhang, T., Shang, E., & Roeder, K.**

Genetic Convergence Analysis of CRISPR Perturbations Deciphers Gene Functional Similarity.

