# locoweed-volatiles-2026

Analysis of volatile organic compound (VOC) profiles in *Oxytropis sericea* (white locoweed) comparing plants with and without the endophytic fungus *Undifilum oxytropis*. Endophyte status is confirmed by seed wash assay (SWA); plants were spatially paired (E+ / E-) to control for microsite and genetic variation.

Data were collected in a greenhouse setting (Spring 2022, Spring 2023) and from natural field populations (2022, 2023). A 2019 historical field dataset is included for supplemental visualization only.

For full methodological detail, analytical decisions, and known data issues, see the [project wiki](../../wiki).

---

## Repository Structure

```
locoweed-volatiles-2026/
├── archive/
│   └── Locoweed All.Rmd
├── data-files/
│   ├── 2022 Field Locoweed.xlsx
│   ├── 2023 Field Locoweed.xlsx
│   ├── Locoweed Field Volatiles 2018.xlsx
│   ├── Locoweed Field Volatiles 2019 Summary.xlsx
│   ├── Spring 2022 GH Locoweed RESULTS.xlsx
│   └── Spring GH 23 Locoweed Volatiles.xlsx
├── figures/
│   ├── Fig1_GH_PCA.png
│   ├── Fig2_Field_PCA.png
│   └── Fig3_Swainsonine.png
├── proc-docs/
│   ├── GH-data-proc.Rmd
│   ├── GH-data-proc.md
│   ├── Field-data-proc.Rmd
│   └── Megan-data-proc.Rmd
├── locoweed_volatiles_analysis.R
├── locoweed-volatiles-2026.Rproj
├── results_summary.txt
└── README.md
```

---

## How to Run

1. Open `locoweed-volatiles-2026.Rproj` in RStudio.
2. Update the absolute data paths in `proc-docs/GH-data-proc.Rmd` to point to your local `data-files/` directory.
3. Run `locoweed_volatiles_analysis.R`. The script calls the three processing Rmds via `knitr::knit()` before executing the main analysis.

---

## Dependencies

Install required packages in R:

```r
install.packages(c("tidyverse", "vegan", "compositions",
                   "ggrepel", "indicspecies", "patchwork"))
```

---

## Authors

Jackson Strand  
Department of Land Resources and Environmental Science  
Montana State University, Bozeman, MT