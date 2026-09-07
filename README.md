# Dianchi hydroacoustic biomass sensitivity

Data and R code supporting the manuscript **“Sensitivity of hydroacoustic fish-biomass estimates to target-strength and length–weight conversions in Dianchi Lake, China”**.

The analysis evaluates how alternative target-strength-to-length (TS–L) equations and length–weight relationships (LWRs) affect hydroacoustic estimates of fish biomass in Dianchi Lake, China. The repository contains the input data and complete R script needed to reproduce the analyses, tables and figures reported in the manuscript.

## Repository structure

```text
.
|-- code/
|   `-- analysis.R
|-- data/
|   |-- hydroacoustic_edsu_data.xlsx
|   |-- fish_catch_data.xlsx
|   |-- species_names.csv
|   `-- table3_dominance.csv
`-- README.md
```

All four files in `data/` are read by `code/analysis.R`. The numerical dominance metrics are recalculated from `fish_catch_data.xlsx`; `table3_dominance.csv` supplies only the species codes, taxonomic families and feeding-guild classifications used in those calculations.

## Software requirements

The verified environment uses R 4.5.1 with these package versions:

- readxl 1.5.0
- ggplot2 4.0.3
- patchwork 1.3.2
- scales 1.4.0
- MCMCglmm 2.36
- ggrepel 0.9.6

Install missing packages before running the analysis:

```r
install.packages(c("readxl", "ggplot2", "patchwork", "scales", "MCMCglmm", "ggrepel"))
```

The analysis script does not install software or alter the user’s package library automatically. It stops with a clear message if a required package is missing.

## How to reproduce the analysis

From the repository root, run:

```bash
Rscript code/analysis.R
```

Alternatively, open the repository folder in RStudio and run:

```r
source("code/analysis.R", encoding = "UTF-8")
```

The script determines the project directory from its own location. Users do not need to modify file paths or add a local `setwd()` call.

The script automatically creates:

- `tables_analysis/`
- `figures_analysis/`
- `session-info.txt`

The two Excel workbooks retain the original Chinese site names, species names and source headers as research records. Their filenames and all executable code are ASCII/English. Immediately after import, the script assigns English analytical variable names and uses `species_names.csv` to map Chinese common names to Latin names.

## Expected verification values

A successful run should report or reproduce:

- 573 valid EDSUs and 14 target-strength classes.
- 2,179 complete individually measured fish used for the length–weight analysis.
- 2,000 moving-block bootstrap replicates with block length 20.
- TS–L sensitivity of approximately 467-fold across the evaluated equations.
- Whole-lake biomass under the Ye et al. (2007) LWR: approximately 36,124 t (95% CI 18,897–61,709 t; CV 32.4%).
- Whole-lake biomass under the community-fitted LWR: approximately 37,805 t (95% CI 19,791–64,107 t; CV 32.3%).
- Community-versus-Ye difference: +4.7% (paired-bootstrap 95% CI 2.6%–7.0%; P = 0.019).
- Size-group correction rates ranging from −36.0% to +36.2%.

Small differences in MCMC summaries can occur if package versions or random-number generators differ. The deterministic tables and seeded bootstrap results should match the verified outputs.

## Data and code availability

The code and data are publicly available in this GitHub repository:

<https://github.com/xiangyy95-hue/dianchi-hydroacoustic-biomass-sensitivity>

A permanent archival version of this project is available from Zenodo:

<https://doi.org/10.5281/zenodo.22476648>

## License and contact

The deposited data and code are available under the [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).

Correspondence may be directed to Chao Guo (`guochao@ihb.ac.cn`) or Wei Li (`liwei@ihb.ac.cn`).
