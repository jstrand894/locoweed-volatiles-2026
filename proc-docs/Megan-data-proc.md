---
title: "Old Volatiles"
author: "Jackson Strand"
date: "2026-04-16"
output: html_document
editor_options: 
  chunk_output_type: console
---


``` r
library(readxl)
library(tidyverse)
```

<!-- READ IN 2018 VOLATILES -->
<!-- OUTPUT: raw_2018 -->

``` r
raw_2018 <- read_excel("/Users/jacksonstrand/Library/CloudStorage/OneDrive-MontanaStateUniversity/Documents/PhD/Locoweed/data-files/Locoweed\ Field\ Volatiles\ 2018.xlsx")
```

```
## Error:
## ! `path` does not exist: '/Users/jacksonstrand/Library/CloudStorage/OneDrive-MontanaStateUniversity/Documents/PhD/Locoweed/data-files/Locoweed Field Volatiles 2018.xlsx'
```

<!-- READ IN 2019 VOLATILES -->
<!-- OUTPUT: raw_2019 -->

``` r
raw_2019 <- read_excel("/Users/jacksonstrand/Library/CloudStorage/OneDrive-MontanaStateUniversity/Documents/PhD/Locoweed/data-files/Locoweed\ Field\ Volatiles\ 2019\ Summary.xlsx", sheet = "ng per g per hr", range = "A1:L41")
```

```
## Error:
## ! `path` does not exist: '/Users/jacksonstrand/Library/CloudStorage/OneDrive-MontanaStateUniversity/Documents/PhD/Locoweed/data-files/Locoweed Field Volatiles 2019 Summary.xlsx'
```


<!-- PROCESS 2019 -->
<!-- OUTPUT: proc_2019 -->

``` r
proc_2019 <-
  raw_2019 %>%
  unite("plant_id", plant, num, sep="_") %>% 
  dplyr::select(-c(biomass)) %>%
  pivot_longer(-c(collection, endophyte, plant_id)) %>%
  mutate(value = ifelse(value == 0, 0.01, value)) %>%
  mutate(value = value + 1) %>%
  pivot_wider(names_from = "name",
              values_from = "value") %>%
    # combine year and plant id into one column
  unite("unite", collection, endophyte, plant_id, sep = "_") %>%
  column_to_rownames("unite") %>%
  clr() %>%
  as.data.frame() %>%
  rownames_to_column("unite") %>%
    # move year and plant id back to being sepearte columns
  separate(unite, 
           into = c("collection", "endophyte", "plant_id", "num"), 
           sep = "_") %>%
  pivot_longer(-c(plant_id, collection, endophyte, num),
               names_to = "name",
               values_to = "value") %>%
  # add a value of 2 to each ngghr to remove negatives
  mutate(value = value + 2) %>%
  pivot_wider(names_from = "name",
              values_from = "value") %>%
  mutate_at(vars(endophyte), list(factor)) %>%
  distinct() %>%
  na.omit()
```

```
## Error:
## ! object 'raw_2019' not found
```



``` r
raw_2018 %>%
  mutate_all(~ replace(., is.na(.), 0))
```

```
## Error:
## ! object 'raw_2018' not found
```

``` r
# oof 2018 has no compound labels just RTs
```













