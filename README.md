# RápidoPGS companion code repository

Last updated: 16/12/2020

This repository contains the code used to run all analyses in the RápidoPGS paper (link to preprint to be added).

### Introduction

For this work, I had to create multiple script versions for testing and debugging, and divide analyses across multiple scripts to be able to run then in the High Performance Computing (HPC) environment. The original scripts I used are in the `original_scripts` directory. This directory contain the most recent versions used in each piece of analysis.
However, to make code easier to read and understand, I created "clean" versions of the code, by grouping related scripts, grouping them, and removing unnecessary content, which can be found in this directory. These scripts should be equivalent to the *original*, though.

### Briefing on each script

Here I'll explain what each script is supposed to do:

* **Creating_RapidoPGS_reference_from_scratch_CLEAN_20201215.R** - This creates and QC a reference panel for RápidoPGS-multi from 1000Genomes Phase III, to be used in the analyses. 
* **Creating_LDpred2_reference_HapMap3_CLEAN_20201215.R** - This creates and QC a reference panel for LDpred2, following LDpred2 author recommendations, including restricting variants to HapMap3 variants only.
* **Preparing_datasets_qcfilt_CLEAN_20201215.R** - This will download and process all 10 datasets used in the paper, following QC guidance recommended by LDpred2 authors, as well as filtering by variants that passed QC in RápidoPGS-multi 1000G panel.
* **Fetch_coordinates.R** - Some datasets (ie. BMI and HEIGHT) come with rsids insead of genomic coordinates. This script will retrieve this information from the panel.
* **rapidoPGS_model_creation_CLEAN_20201215.R** - This runs all RápidoPGS analyses, including RápidoPGS-single and RápidoPGS-multi (auto and informed prior, and both sets of alpha parameters) for all traits.
* **LDpred2_model_creation_CLEAN_20201215.R** - This runs all LDpred2 analyses. Due to time constraints, we made this script run using a 1-10 number argument corresponding to each dataset, so each dataset can be run individually. 
* **Compute_custom_priors_BMI_Height_CLEAN_20201215.R** - This computes the informed prior for both BMI and Height datasets, used in RápidoPGS analyses.
* **Benchmarking_rapidoPGS_LDpred2_CLEAN_20201215.R** - This benchmarks LDpred2 and RápidoPGS runtime for the Asthma dataset (corresponding to Figure 2 in the paper).
* **Benchmarking_rapidoPGS_alltraits_CLEAN_20201215.R** - This benchmarks RápidoPGS flavours runtime for all datasets (corresponding to Figure S1 in the paper).


*Guillermo Reales*
