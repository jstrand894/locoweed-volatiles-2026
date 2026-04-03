---
title: "Jackson Locoweed GH Data Processing"
author: "Jackson Strand"
date: "2026-03-17"
output: html_document
editor_options: 
  chunk_output_type: console
---


``` r
library(compositions)
library(tidyverse)
library(readxl)
```

<!-- READ IN RAW SPRING 2023 DATA -->
<!-- OUTPUT: gh22spring & gh23spring -->

``` r
# SPRING 2022 DATA
gh22spring <-
read_excel("/Users/jacksonstrand/Library/CloudStorage/OneDrive-MontanaStateUniversity/MSU\ PhD/Projects/Locoweed/data-files/Spring\ 2022\ GH\ Locoweed\ RESULTS.xlsx", sheet = "RAW2")
```

```
## Error:
## ! `path` does not exist: '/Users/jacksonstrand/Library/CloudStorage/OneDrive-MontanaStateUniversity/MSU PhD/Projects/Locoweed/data-files/Spring 2022 GH Locoweed RESULTS.xlsx'
```

``` r
# SPRING 2023 DATA
gh23spring <-
  read_excel("/Users/jacksonstrand/Library/CloudStorage/OneDrive-MontanaStateUniversity/MSU\ PhD/Projects/Locoweed/data-files/Spring\ GH\ 23\ Locoweed\ Volatiles.xlsx", sheet = "RAW_PROC")
```

```
## Error:
## ! `path` does not exist: '/Users/jacksonstrand/Library/CloudStorage/OneDrive-MontanaStateUniversity/MSU PhD/Projects/Locoweed/data-files/Spring GH 23 Locoweed Volatiles.xlsx'
```

<!-- CLR ALL SPRING DATA -->
<!-- OUTPUT: joined.data.clr  -->

``` r
joined.data <-
gh23spring %>%
  dplyr::select(plant.id, compound, ngghr) %>%
  # filter(compound
  group_by(plant.id, compound) %>%
  na.omit() %>%
  dplyr::summarise(plant.id, compound, ngghr = sum(ngghr)) %>%
  distinct() %>%
  pivot_wider(names_from = "compound",
              values_from = "ngghr") %>%
  pivot_longer(-plant.id,
               names_to = "compound",
               values_to = "ngghr") %>%
  mutate(year = "Spring 2023") %>%
  left_join(gh23spring %>%
              dplyr::select(plant.id, 
                            swa, pair),
            by = "plant.id") %>%
  distinct() %>%
  rbind(gh22spring %>%
          dplyr::select(plant.id, compound, ngghr) %>%
          # filter(compound
          group_by(plant.id, compound) %>%
          na.omit() %>%
          dplyr::summarise(plant.id, compound, ngghr = sum(ngghr)) %>%
          distinct() %>%
          pivot_wider(names_from = "compound",
                      values_from = "ngghr") %>%
          pivot_longer(-plant.id,
                       names_to = "compound",
                       values_to = "ngghr") %>%
          distinct() %>%
          mutate(year = "Spring 2022") %>%
          left_join(gh22spring %>%
                      dplyr::select(plant.id, 
                                    swa, pair),
                    by = "plant.id")) %>%
        mutate_at(vars(pair, swa), list(factor)) %>%
  distinct()
```

```
## Error in `dplyr::summarise()`:
## ℹ In argument: `plant.id`.
## ℹ In group 3: `plant.id = "10A"`, `compound =
##   "decanal"`.
## Caused by error:
## ! `plant.id` must be size 1, not 2.
## ℹ To return more or less than 1 row per group,
##   use `reframe()`.
```

``` r
# joined.data %>%
#   filter(plant.id == "10A", compound == "decanal")

gh.joined.data.clr <-
  joined.data %>%
  dplyr::select(plant.id, ngghr, compound, year) %>%
  mutate(ngghr = ifelse(ngghr == 0, 0.01, ngghr)) %>%
  mutate(ngghr = ngghr + 1) %>%
  pivot_wider(names_from = "compound",
              values_from = "ngghr") %>%
    # combine year and plant id into one column
  unite(plant.id, plant.id, year, sep = "_") %>%
  column_to_rownames("plant.id") %>%
  clr() %>%
  as.data.frame() %>%
  rownames_to_column("plant.id") %>%
    # move year and plant id back to being sepearte columns
  separate(plant.id, into = c("plant.id", "year"), sep = "_") %>%
  pivot_longer(-c(plant.id, year),
               names_to = "compound",
               values_to = "ngghr") %>%
  # add a value of 2 to each ngghr to remove negatives
  mutate(ngghr = ngghr + 2) %>%
  pivot_wider(names_from = "compound",
              values_from = "ngghr") %>%
  left_join(joined.data %>%
              dplyr::select(plant.id, swa, year, pair) %>%
              mutate(plant.id = as.character(plant.id)),
            by = c("plant.id","year")) %>%
  mutate_at(vars(swa, year, pair), list(factor)) %>%
  distinct() %>%
  na.omit()
```






