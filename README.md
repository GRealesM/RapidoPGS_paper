## Structure of the code


#### Prepare data

- `Preparing_datasets_qcfilt_clean_20210317.R`: Download datasets, apply QC, generate QC plots (Figures S1-10), filter them by the HapMap3 variants, and save them. 

- `Creating_LDpred2_reference_HapMap3_20201111.R`: Download UKBB pre-computed LD matrices of European individuals, to be used with LDpred2 and RápidoPGS.

- `prscs_preprocess_20210317.R`: Preprocess summary statistics to satisfy PRScs format and requirements.

- `sbayesr_preprocess_20210317.R`: Pre-process summary statistics to satisfy SBayesR format and requirements.


#### Run methods

- `rapidopgs_run_20210317.R`: Run RápidoPGS (R script, to be called from slurm).
- `LDpred2_gwide_run_20210223.R`: Run LDpred2-auto genome-wide (R script, to be called from slurm). 
- `LDpred2_perchr_run_20210317.R`: Run LDpred2-auto genome-wide (R script, to be called from slurm). LDpred2 results in the paper correspond to the per-chromosome approach.
- `slurm_run_ldpred2perchr`: Slurm script for LDpred2-auto per-chromosome approach.
- `slurm_run_PRScs`: Run PRScs (slurm script).
- `slurm_run_rapidopgs: Run RápidoPGS (slurm script to run R code on the HPC).
- `slurm_run_SBayesR`: Run SBayesR (slurm script).


#### Timings

- `Benchmarking_all_20210317.R`: Time all methods except for PRScs, which is timed independently (R script, to be called from slurm).
- `slurm_run_mbm_all`: Slurm script for timing all methods (except for PRScs).
- `slurm_run_mbm_PRScs`: Slurm script for timing PRScs.

**NOTE**: PRScs runs in Python and couldn't be called from R, so it was run (and timed) straight from slurm. See relevant code for details.


#### Visualise results

- `Figures_workshop.R`: Create figures for manuscript (Figures 1-4) and compute differences and ratios for models AUC and r2.


#### Helper scripts

- `Compute_custom_priors_BMI_Height_20210223.R`: Compute informed priors to be used for RápidoPGS-multi approach.
- `Compute_freqs.sh`: Compute allele frequencies from CEU population in 1000 Genomes, if they're not provided with the dataset.
- `Mergeback.R`: Merge everything after computing frequencies.
- `Fetch_coordinates.R`: Use 1000 Genomes to get hg19 coordinates from SNPs in the summary statistics dataset, if missing.


