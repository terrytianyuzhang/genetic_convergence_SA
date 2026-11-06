# Persistent experiment decisions

- Use Xin `151 TH Prkcd Grin2c Glut` cells with `gene_target == "Non_target"`.
- Use only the supplied 10-celltype/10-perturbation h5ad; no dataset concatenation.
- Rank genes by total raw expression across ALL eligible 151 controls before sampling cells.
- Keep the top 1500 genes; ties follow original gene-column order. Consecutive blocks of 50 genes form 30 modules.
- Sort eligible cell IDs deterministically, sample 1600 unique cells with seed 4, allocate A-F 200 each and NT 400.
- Keep the original A-F module multipliers and original 80% gene dilution, with mask seeds A=1 through F=6.
- Multiply the selected raw counts by the gene-specific multiplier and then use R `round()`. NT remains unchanged.
- **Always draw union AND intersect. For the paper, prefer INTERSECT.**
- **The convergence pair p-value is MAX of the forward and reverse p-values.**
- Keep the original convergence_testing parameters: dense_pca, group_lasso, lambda.1se, 10 folds, standardize_feature=TRUE.
- Keep the original Poisson GLM and pooled-standardized Chen benchmark.
- Each HMC pair/direction has a fixed separate seed, 10001-10030. Each Chen group/module also has a fixed seed, 50001-50180.
- All later stages load the SAME `data/Xin_151_200_400.Rdata` from script 01. They never resample cells or genes.
- Recomputed data and numerical test results are compared to existing saved objects and stop if different. Input checksums use base R tools::md5sum().
- Exact numerical reproducibility assumes the same server R/package/numerical-library environment. Each script records sessionInfo; no extra environment-lock or thread-controller package is required.
- **No local HMC/Chen computations.** Scripts 02 and 04 run on the server. Script 01 now defaults to LOCAL paths at the user's explicit request; preserve its server input/output paths as comments. Editing paths does not authorize the assistant to run the experiment locally.
- Both the local and server project directories are named **redo_by_Xin_data**, with Xin, not XIn. Do not attribute file-not-found errors to an unverified capitalization difference.
- **All five scripts should resemble the old folder:** library() calls, explicit paths/parameters, sequential code. Avoid cfg frameworks, xin_main entry points, automatic R restarts, mandatory thread-control packages, and additional checksum packages.

# Files and packages

Script 01 active local paths:
- Input: `/Users/erganshang/Desktop/pert_map/Xin/data/raw/subset_adata_SCARF_10celltypes_10perturbations.h5ad`
- Data: `/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/data/`
- Figures: `/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/figures/`

Server project root:
`/raid6/home/eshang/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data`

- Input: `raw/subset_adata_SCARF_10celltypes_10perturbations.h5ad`
- Shared data: `01/data/Xin_151_200_400.Rdata`
- HMC results: `01/data/convergence/`
- Chen results: `01/data/Chen/`
- Figures: `01/figures/`

Packages already used by the original folder:
- 01: ggplot2, reshape2, withr (the old helper already used withr::with_seed).
- 02: HMC.
- 03: HMC, ggplot2.
- 04: highmean.
- 05: ggplot2.

The only new direct package is **hdf5r**, needed only in 01 to read h5ad. There is no direct digest or RhpcBLASctl dependency.

# Running

Script 01 defaults to local input and local output. Its three server path assignments remain commented out for switching back. If running 01 locally, copy the resulting `Xin_151_200_400.Rdata` unchanged into server `01/data/` before running 02-05. Alternatively, activate the three commented server paths in 01 and run the full workflow on the server. Scripts 02-05 retain their server paths.
Restart the current R session after the previous recursive errors, then run in order 01 -> 05 in the chosen environments.
`source()` now executes each script directly, like the original folder; do not call the removed xin_main_* functions.
Alternatively use Rscript for each file, in order. Do not run concurrent copies of a stage against the same output directory.

The updated scripts use `Xin_151_200_400.Rdata`, not the early framework's `prepared_data.rds`. Run revised 01 first and use revised 02-05 together. Do not mix the early framework scripts/manifests/results with this version. Existing old outputs are not deleted automatically.

Data output files are deterministic under the same environment. PDF visual content is deterministic; standard ggsave PDF creation metadata may differ between runs, as in the original scripts.

Multiplication/rounding preserves the effect design. Actual observed treatment-vs-NT fold changes and power can differ from normal simulations; zero counts remain zero.
