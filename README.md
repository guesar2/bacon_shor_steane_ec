# README

This repository contains the source code used to generate the data for the paper _Improved Performance of the Bacon-Shor Code with Steane's Syndrome Extraction Method_ by Guillermo Escobar and Mauricio Gutiérrez ([https://arxiv.org/abs/2403.01659](https://arxiv.org/abs/2403.01659)).

## Overview
The code uses a simulation toolkit with [CHP](https://arxiv.org/abs/quant-ph/0406196) as its core simulator. The main simulation scripts are located in the `src` directory:
- `sim_shor_ec.py`
- `sim_steane_ec.py`

These scripts generate JSON files containing the results of the sampling. Since subset sampling is employed, a separate JSON file is generated for each subset.

## Calculating Logical Error Rates
To compute logical error rates:
1. Summarize the generated JSON files into a single JSON file using the `summarize.py` script.
2. Use the `plogical.py` script to calculate logical error rates from the summarized data.

## Automation with Bash Scripts
The following bash scripts automate the simulation and result generation process:
- `run_simulations_shor.sh`
- `run_simulations_steane.sh`
- `generate_results_shor.sh`
- `generate_results_steane.sh`

Running these scripts should generate most of the data used in the paper. However, since the simulations may take several days to complete, it is recommended to adjust the variables and run smaller, more manageable parameter sets iteratively. Also, note that the last two scripts require the jq JSON processor.
