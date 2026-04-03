# ============================================================
# locoweed_volatiles_analysis.R
# Locoweed VOC Profile Analysis: Endophyte Effect
# Author: Jackson Strand
# ============================================================
#
# RESEARCH QUESTION:
#   Do locoweed plants (Oxytropis sericea) harboring the endophytic
#   fungus Undifilum oxytropis (SWA-positive / E+) differ in their
#   volatile organic compound (VOC) profiles from plants without the
#   endophyte (SWA-negative / E-)? Plants were spatially paired
#   (one E+, one E-) to control for microsite and genetic variation.
#
# DATASETS:
#   - Greenhouse (Jackson): Spring 2022 and Spring 2023
#     Whole-plant VOC collection; endophyte status confirmed by SWA
#     assay; plants paired by seed family / planting block.
#   - Field (Jackson): 2022 and 2023
#     Flowers and leaves sampled separately from co-occurring paired
#     plants at a natural population.
#   - Historical field (Megan): 2019
#     Used for supplemental visualization only; not included in
#     main PERMANOVA (different compound set, no pairing metadata).
#
# STATISTICAL APPROACH:
#   1. CLR transformation (centered log-ratio) on raw VOC ng/g/hr
#      values, applied once to pseudocount-adjusted data.
#      Euclidean distance on CLR data = Aitchison distance, which
#      is the appropriate metric for compositional VOC profiles.
#   2. PERMANOVA (adonis2, vegan) to test multivariate VOC differences.
#      - Terms: swa + year + pair (greenhouse) or
#               swa + collection.type + year + pair2 (field)
#      - by = "margin": each term tested after all others (Type III),
#        appropriate because we want the effect of endophyte AFTER
#        accounting for year and pairing.
#      - strata = year: permutations restricted within years to
#        account for the non-exchangeability of observations across years.
#   3. Betadisper: test homogeneity of within-group dispersion.
#      Significant dispersion differences can inflate PERMANOVA p-values.
#   4. PCA for visualization of compositional space.
#   5. SIMPER: which compounds most contribute to group dissimilarity.
#   6. IndVal: compounds significantly associated with E+ or E-.
#
# ============================================================
# AUDIT NOTES — ISSUES FOUND IN Locoweed All.Rmd
# ============================================================
#
# CRITICAL:
#
#   1. DOUBLE (AND TRIPLE) CLR IN MAIN PERMANOVA CHUNKS
#      In Field-data-proc.Rmd, field_comb_clean applies log2()
#      BEFORE CLR, so field_comb_clr_final contains
#      CLR(log2(x)). Then in the Locoweed All.Rmd PERMANOVA
#      chunks, the data is further transformed as (ngghr^2)*100000
#      and CLR() is applied again — resulting in a triple
#      transformation with no compositional interpretation.
#      GH-data-proc.Rmd correctly applies CLR to raw values, but
#      the main analysis then squares and re-CLRs those too.
#      FIX: This script applies CLR exactly once, to raw values
#      from joined.data and field_comb.
#
#   2. WRONG PCA OBJECT IN GH AXIS LABELS
#      The GH PCA plot references summary(pca_result) (the FIELD
#      pca object) for axis label percentages instead of
#      summary(pca_result_gh). The GH plot shows the wrong
#      variance explained values.
#      FIX: Script uses the correct GH pca object for GH labels.
#
#   3. pca_data vs pca_data_field MISMATCH
#      The "FIELD PERMANOVA + PCA" chunk computes an object called
#      pca_data, but the FIELD_PCA_Plot ggplot call uses
#      pca_data_field. The plot relies on a leftover object from
#      a previous session run — running the chunk fresh would fail
#      or silently plot the wrong data.
#      FIX: This script creates and uses consistently named objects.
#
#   4. UNDOCUMENTED MANUAL CORRECTIONS IN BAR PLOT DATA
#      The "field volatiles data" chunk applies unexplained numeric
#      correction factors to 2022 data:
#        - 2022 leaves: mean_ngghr/3 and se/3
#        - 2022 flowers: mean_ngghr/22 and se/22
#        - 2022 bocimene (flowers): mean_ngghr*19 and se*5
#      No comments explain what these are correcting for
#      (e.g., unit conversions, collection volume differences).
#      These should be traced back to a documented source and
#      applied at the data processing stage, not silently in a
#      plot pipeline.
#      ACTION NEEDED: Verify what these correction factors represent
#      before publication.
#
#   5. dlimonene RELABELED AS linalool IN FIELD PLOT
#      In the field volatiles plot chunk:
#        mutate(compound = if_else(compound == "dlimonene", "linalool", compound))
#      These are chemically distinct compounds. This relabeling
#      is either a compound ID correction or an error.
#      ACTION NEEDED: Confirm whether the compound identified as
#      dlimonene in the field data is actually linalool.
#
#   6. HARDCODED SUMMARY VALUES ADDED VIA bind_rows
#      Several rows are manually inserted into field_data with
#      hardcoded mean_ngghr and se values (e.g., dlimonene 2023
#      E-, mean=6.21, se=1.5). No source or derivation is given.
#      These values should be derived programmatically from the
#      raw data.
#      ACTION NEEDED: Replace hardcoded values with calculations.
#
# MODERATE:
#
#   7. INCONSISTENT by ARGUMENT IN PERMANOVA
#      GH PERMANOVA uses by = "terms" (sequential / Type I).
#      Field PERMANOVA uses by = "margin" (marginal / Type III).
#      For the question "does endophyte status matter after
#      controlling for year and pair?", by = "margin" is correct
#      for both. Sequential results depend on term order and do
#      not cleanly isolate the endophyte effect.
#      FIX: This script uses by = "margin" for both.
#
#   8. strata = year VS strata = pair
#      For a paired design, restricting permutations within pairs
#      (strata = pair) would be most appropriate — it directly
#      preserves the pairing structure. strata = year accounts
#      for year as a block but does not enforce pairing.
#      This is a minor issue given sample sizes; current approach
#      is a reasonable conservative choice and is kept here.
#
#   9. OLD DUPLICATE PERMANOVA CODE STILL IN Locoweed All.Rmd
#      An older "PERMANOVA_FIELD" chunk precedes the clean
#      "FIELD PERMANOVA + PCA" chunk and uses for_permanova
#      (wrong object) for building dist.matrix.field. Running
#      both chunks would produce conflicting results.
#
#   10. PATH HARDCODING IN GH-data-proc.Rmd
#       GH-data-proc.Rmd reads from absolute paths pointing to a
#       different OneDrive directory than the current project root.
#       This will break if run on any other machine or if the
#       directory is reorganized. Field-data-proc.Rmd correctly
#       uses relative paths.
#
# MINOR:
#
#   11. INCOMPLETE CHUNK
#       The "Greenhouse Table" chunk at the top of Locoweed All.Rmd
#       contains only "joined.data %>%" with no further code.
#
#   12. INCONSISTENT VARIABLE NAMING
#       Mix of dot-separated (permanova.output.field) and
#       underscore-separated (permanova_output_gh) names throughout.
#
#   13. set.seed(1)
#       Seed value of 1 is fine functionally but using a more
#       arbitrary value reduces the appearance of seed-searching.
#
# ============================================================

library(tidyverse)
library(vegan)
library(compositions)
library(ggrepel)
library(indicspecies)
library(patchwork)


# ============================================================
# SECTION 1: LOAD DATA VIA PROCESSING SCRIPTS
# ============================================================
# The processing Rmd files handle data import and light cleaning.
# Key outputs produced:
#   joined.data        - Long-format GH VOC data, raw ng/g/hr (both years)
#   field_comb         - Long-format field VOC data, raw ng/g/hr (both years)
#   proc_2019          - CLR-processed 2019 field data (Megan; visualization only)
#
# NOTE: GH-data-proc.Rmd also produces gh.joined.data.clr and
# field_comb_clr_final, but this script re-derives CLR correctly
# from the raw objects above, so those outputs are not used here.
#
# NOTE: GH-data-proc.Rmd uses absolute paths to a separate OneDrive
# directory. Update those paths to point to data-files/ in this
# project if running on a new machine or after reorganizing folders.

knitr::knit("proc-docs/GH-data-proc.Rmd",    envir = globalenv())
knitr::knit("proc-docs/Megan-data-proc.Rmd", envir = globalenv())
knitr::knit("proc-docs/Field-data-proc.Rmd", envir = globalenv())


# ============================================================
# SECTION 2: GREENHOUSE ANALYSIS — SPRING 2022 & 2023
# ============================================================

# Compounds to exclude: internal standards, known contaminants, and non-biogenic compounds.
# Synthetic/industrial: homosalate, ethylhexylsalicylate, diethylphthalate (UV filters /
#   plasticizers), benzothiazole (rubber additive), butylatedhydroxtoluene / BHT (synthetic
#   antioxidant; common SPME fiber contaminant), dodecylacrylate (synthetic acrylate),
#   caprylyl acetate (not established as biogenic in Oxytropis).
# Straight-chain alkane contaminants: hexadecane, tridecane, pentadecane, nanodecane.
# Anthropogenic aromatic: p-trimethyl (likely trimethylbenzene).
gh_exclude <- c("tetradecane", "tetradec", "tetramethyl",
                "standard", "o+C5xe", "cyclohexene",
                "homosalate", "ethylhexylsalicylate", "diethylphthalate",
                "benzothiazole", "butylatedhydroxtoluene",
                "dodecylacrylate", "caprylyl acetate",
                "hexadecane", "tridecane", "pentadecane", "nanodecane",
                "p-trimethyl", "dimethylnona")

# --- 2A. Build compound matrix from raw data ---
# Aggregate multiple measurements per plant × compound to one value.
gh_wide <- joined.data %>%
  filter(!compound %in% gh_exclude,
         !is.na(ngghr)) %>%
  group_by(plant.id, year, swa, pair, compound) %>%
  summarise(ngghr = sum(ngghr), .groups = "drop") %>%
  pivot_wider(names_from  = compound,
              values_from = ngghr,
              values_fill = 0) %>%
  arrange(plant.id, year)

# Separate metadata from compound abundances
gh_meta <- gh_wide %>%
  select(plant.id, year, swa, pair) %>%
  mutate(row_id = paste(plant.id, year, pair, sep = "~")) %>%
  column_to_rownames("row_id")

gh_comp_mat <- gh_wide %>%
  select(-plant.id, -year, -swa, -pair) %>%
  as.matrix()
rownames(gh_comp_mat) <- rownames(gh_meta)

# Replace true zeros with a small pseudocount before CLR.
# CLR requires strictly positive values; 0.01 is consistent with
# the original processing pipeline.
gh_comp_mat[gh_comp_mat == 0] <- 0.01

# --- 2B. CLR transformation (applied once, to raw values) ---
# CLR(x)_i = log(x_i / geometric_mean(x))
# Euclidean distance on CLR data = Aitchison distance.
gh_clr <- clr(gh_comp_mat)

# Drop all-zero rows and samples missing swa/pair (e.g., plant 3A2 Spring 2023)
nonzero_gh <- rowSums(as.matrix(gh_clr)) != 0
gh_meta    <- gh_meta |> filter(nonzero_gh, !is.na(swa), !is.na(pair))
gh_clr     <- gh_clr[rownames(gh_meta), ]

# --- 2C. Aitchison distance matrix ---
gh_dist <- vegdist(as.matrix(gh_clr), method = "euclidean")

# --- 2D. Betadisper: test homogeneity of dispersion ---
# PERMANOVA p-values can be inflated if within-group spread differs
# between groups (i.e., the groups differ in variance, not just
# centroid location). Test for year and swa separately.
bd_gh_year <- betadisper(gh_dist, gh_meta$year)
bd_gh_swa  <- betadisper(gh_dist, gh_meta$swa)

get_betadisp_p <- \(bd) permutest(bd)$tab |>
  as.data.frame() |>
  rownames_to_column("term") |>
  filter(term == "Groups") |>
  pull(`Pr(>F)`)

betadisp_gh_year_p <- get_betadisp_p(bd_gh_year)
betadisp_gh_swa_p  <- get_betadisp_p(bd_gh_swa)

cat("GH betadisper by year, p =", betadisp_gh_year_p, "\n")
cat("GH betadisper by swa,  p =", betadisp_gh_swa_p,  "\n")
# If either p < 0.05, interpret the PERMANOVA result cautiously —
# the significant term may partly reflect dispersion, not location.

# --- 2E. PERMANOVA ---
# Tests whether VOC profiles differ by endophyte status (swa),
# after accounting for plant pair identity.
# by = "margin": marginal (Type III) test — order of terms does not
#   affect results.
# strata = year: permutations restricted within years to control for
#   year as a blocking factor. Year is NOT included as a formula term
#   because it is fully collinear with strata and cannot be estimated
#   as a separate term (produces Df = 0).
set.seed(4721)
gh_permanova <- adonis2(
  gh_dist ~ swa + pair,
  data         = gh_meta,
  permutations = 999,
  strata       = gh_meta$year,
  by           = "margin"
)
print(gh_permanova)

# --- 2E-ii. Year-stratified PERMANOVAs ---
# Run separate models per year to confirm consistency of swa effect
# across years and to allow pair to be tested within-year.
run_gh_year_permanova <- function(year_val) {
  keep <- gh_meta$year == year_val
  d    <- vegdist(as.matrix(gh_clr)[keep, ], method = "euclidean")
  set.seed(4721)
  adonis2(d ~ swa + pair, data = gh_meta[keep, ],
          permutations = 999, by = "margin")
}

gh_perm_2022 <- run_gh_year_permanova("Spring 2022")
gh_perm_2023 <- run_gh_year_permanova("Spring 2023")

cat("GH PERMANOVA Spring 2022:\n"); print(gh_perm_2022)
cat("GH PERMANOVA Spring 2023:\n"); print(gh_perm_2023)

# --- 2F. PCA for visualization ---
gh_pca    <- prcomp(as.matrix(gh_clr), scale. = TRUE)
gh_scores <- as.data.frame(gh_pca$x) %>%
  rownames_to_column("row_id")

gh_pca_data <- gh_meta %>%
  rownames_to_column("row_id") %>%
  left_join(gh_scores, by = "row_id") %>%
  mutate(
    swa  = factor(swa,  levels = c(0, 1), labels = c("Absent", "Present")),
    year = recode(year, "Spring 2022" = "2022", "Spring 2023" = "2023")
  )

# Variance explained by PC1 and PC2
gh_pca_imp <- summary(gh_pca)$importance |>
  as.data.frame() |>
  slice(2)  # Proportion of Variance row
gh_pc1_var <- round(gh_pca_imp$PC1 * 100, 1)
gh_pc2_var <- round(gh_pca_imp$PC2 * 100, 1)

# Top loading compounds (useful for interpretation)
gh_rot <- as.data.frame(gh_pca$rotation) |> rownames_to_column("compound")
gh_top_pc1 <- gh_rot |> slice_max(abs(PC1), n = 3) |> pull(compound)
gh_top_pc2 <- gh_rot |> slice_max(abs(PC2), n = 3) |> pull(compound)
cat("GH PC1 top compounds:", paste(gh_top_pc1, collapse = ", "), "\n")
cat("GH PC2 top compounds:", paste(gh_top_pc2, collapse = ", "), "\n")

# Loading arrows: scale loadings for biplot overlay
gh_pca_loadings <- as.data.frame(gh_pca$rotation[, 1:2]) %>%
  rownames_to_column("compound") %>%
  mutate(
    PC1 = PC1 * max(abs(gh_pca_data$PC1)) * 0.8,
    PC2 = PC2 * max(abs(gh_pca_data$PC2)) * 0.8
  ) %>%
  # Keep only the most influential compounds for a clean biplot
  slice_max(order_by = sqrt(PC1^2 + PC2^2), n = 6)

GH_PCA_Plot <-
  ggplot(gh_pca_data) +
  geom_point(aes(PC1, PC2, shape = swa, fill = year),
             size = 2, stroke = 0.3, alpha = 0.7, color = "black") +
  stat_ellipse(aes(PC1, PC2, color = year, group = year), level = 0.95) +
  geom_segment(data = gh_pca_loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red") +
  geom_text_repel(data = gh_pca_loadings,
                  aes(PC1, PC2, label = compound),
                  color = "black", size = 3,
                  box.padding = 0.4, max.overlaps = Inf) +
  scale_fill_manual(values  = c("2022" = "#402747", "2023" = "#BFD8B8")) +
  scale_color_manual(values = c("2022" = "#402747", "2023" = "#BFD8B8")) +
  scale_shape_manual(values = c("Absent" = 21, "Present" = 25)) +
  labs(
    x     = paste0("PC1 (", gh_pc1_var, "%)"),
    y     = paste0("PC2 (", gh_pc2_var, "%)"),
    fill  = "Year",
    color = "Year",
    shape = "Endophyte"
  ) +
  theme_bw() +
  theme(
    panel.grid       = element_blank(),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  )

GH_PCA_Plot


# ============================================================
# SECTION 3: FIELD ANALYSIS — 2022 & 2023
# ============================================================
# Core compounds detected consistently across both years.
# This pre-selection focuses the analysis on compounds that are
# biologically meaningful and reliably detected; rare or trace
# compounds are excluded to reduce noise in the distance matrix.
field_compounds <- c("heptanal", "sixmethyl", "octanal", "z3hex",
                     "dlimonene", "nonanal", "decanal",
                     "methylsalicylate", "cisthreehexiso", "bocimene")

# --- 3A. Build compound matrix from raw data ---
# Exclude non-tissue samples (type "s" = senesced, 0 = unknown)
field_wide <- field_comb %>%
  filter(
    compound %in% field_compounds,
    !is.na(collection.type),
    !collection.type %in% c("senesced_leaves", "0")
  ) %>%
  group_by(plant.id, year, swa, collection.type, pair2, compound) %>%
  summarise(ngghr = sum(ngghr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = compound,
              values_from = ngghr,
              values_fill = 0)

field_meta <- field_wide %>%
  select(plant.id, year, swa, collection.type, pair2) %>%
  mutate(row_id = paste(plant.id, year, collection.type, pair2, sep = "~")) %>%
  column_to_rownames("row_id")

field_comp_mat <- field_wide %>%
  select(all_of(field_compounds)) %>%
  as.matrix()
rownames(field_comp_mat) <- rownames(field_meta)

field_comp_mat[field_comp_mat == 0] <- 0.01

# --- 3B. CLR transformation ---
field_clr     <- clr(field_comp_mat)
nonzero_field <- rowSums(as.matrix(field_clr)) != 0
field_meta    <- field_meta |> filter(nonzero_field)
field_clr     <- field_clr[rownames(field_meta), ]

# --- 3C. Aitchison distance matrix ---
field_dist <- vegdist(as.matrix(field_clr), method = "euclidean")

# --- 3D. Betadisper ---
bd_field_year <- betadisper(field_dist, field_meta$year)
bd_field_type <- betadisper(field_dist, field_meta$collection.type)
bd_field_swa  <- betadisper(field_dist, field_meta$swa)

betadisp_field_year_p <- get_betadisp_p(bd_field_year)
betadisp_field_type_p <- get_betadisp_p(bd_field_type)
betadisp_field_swa_p  <- get_betadisp_p(bd_field_swa)

cat("Field betadisper by year,            p =", betadisp_field_year_p, "\n")
cat("Field betadisper by collection type, p =", betadisp_field_type_p, "\n")
cat("Field betadisper by swa,             p =", betadisp_field_swa_p,  "\n")

# --- 3E. PERMANOVA ---
# collection.type (flowers vs leaves) is included to ask whether
# endophyte status matters after accounting for tissue differences.
# Year dropped from formula — collinear with strata (see GH note above).
set.seed(4721)
field_permanova <- adonis2(
  field_dist ~ swa + collection.type + pair2,
  data         = field_meta,
  permutations = 999,
  strata       = field_meta$year,
  by           = "margin"
)
print(field_permanova)

# --- 3E-ii. Year-stratified PERMANOVAs ---
run_field_year_permanova <- function(year_val) {
  keep <- field_meta$year == year_val
  d    <- vegdist(as.matrix(field_clr)[keep, ], method = "euclidean")
  set.seed(4721)
  adonis2(d ~ swa + collection.type + pair2, data = field_meta[keep, ],
          permutations = 999, by = "margin")
}

field_perm_2022 <- run_field_year_permanova("2022")
field_perm_2023 <- run_field_year_permanova("2023")

cat("Field PERMANOVA 2022:\n"); print(field_perm_2022)
cat("Field PERMANOVA 2023:\n"); print(field_perm_2023)

# --- 3F. PCA for visualization ---
# Combined PCA, faceted by tissue type in the plot.
field_pca    <- prcomp(as.matrix(field_clr), scale. = TRUE)
field_scores <- as.data.frame(field_pca$x) %>%
  rownames_to_column("row_id")

field_pca_data <- field_meta %>%
  rownames_to_column("row_id") %>%
  left_join(field_scores, by = "row_id") %>%
  mutate(
    swa             = factor(swa, levels = c(0, 1), labels = c("Absent", "Present")),
    collection.type = recode(collection.type, "flowers" = "Flowers", "leaves" = "Leaves")
  )

field_pca_imp <- summary(field_pca)$importance |>
  as.data.frame() |>
  slice(2)  # Proportion of Variance row
field_pc1_var <- round(field_pca_imp$PC1 * 100, 1)
field_pc2_var <- round(field_pca_imp$PC2 * 100, 1)

field_rot <- as.data.frame(field_pca$rotation) |> rownames_to_column("compound")
field_top_pc1 <- field_rot |> slice_max(abs(PC1), n = 3) |> pull(compound)
field_top_pc2 <- field_rot |> slice_max(abs(PC2), n = 3) |> pull(compound)
cat("Field PC1 top compounds:", paste(field_top_pc1, collapse = ", "), "\n")
cat("Field PC2 top compounds:", paste(field_top_pc2, collapse = ", "), "\n")

field_pca_loadings <- as.data.frame(field_pca$rotation[, 1:2]) %>%
  rownames_to_column("compound") %>%
  mutate(
    PC1 = PC1 * max(abs(field_pca_data$PC1)) * 0.8,
    PC2 = PC2 * max(abs(field_pca_data$PC2)) * 0.8
  ) %>%
  slice_max(order_by = sqrt(PC1^2 + PC2^2), n = 6)

FIELD_PCA_Plot <-
  ggplot(field_pca_data) +
  geom_point(aes(PC1, PC2, shape = swa, fill = year),
             size = 2, stroke = 0.3, alpha = 0.8, color = "black") +
  stat_ellipse(aes(PC1, PC2, color = year, group = year), level = 0.95) +
  geom_segment(data = field_pca_loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red") +
  geom_text_repel(data = field_pca_loadings,
                  aes(PC1, PC2, label = compound),
                  color = "black", size = 3,
                  box.padding = 0.4, max.overlaps = Inf) +
  scale_fill_manual(values  = c("2022" = "#402747", "2023" = "#BFD8B8")) +
  scale_color_manual(values = c("2022" = "#402747", "2023" = "#BFD8B8")) +
  scale_shape_manual(values = c("Absent" = 21, "Present" = 25)) +
  facet_wrap(~ collection.type) +
  labs(
    x     = paste0("PC1 (", field_pc1_var, "%)"),
    y     = paste0("PC2 (", field_pc2_var, "%)"),
    fill  = "Year",
    color = "Year",
    shape = "Endophyte"
  ) +
  theme_bw() +
  theme(
    panel.grid        = element_blank(),
    strip.background  = element_rect(fill = NA),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  )

FIELD_PCA_Plot


# ============================================================
# SECTION 4: SIMPER — COMPOUND CONTRIBUTIONS TO DISSIMILARITY
# ============================================================
# SIMPER asks: which compounds most explain the dissimilarity
# between endophyte-present and endophyte-absent groups?

# --- 4A. SIMPER on raw ng/g/hr values (Bray-Curtis) ---
# CLR-transformed values can be negative, which breaks Bray-Curtis.
# Raw ng/g/hr values are non-negative and appropriate for simper().
# Build raw matrices aligned with the already-filtered metadata.

gh_raw_mat <- gh_wide |>
  mutate(row_id = paste(plant.id, year, pair, sep = "~")) |>
  filter(row_id %in% rownames(gh_meta)) |>
  arrange(match(row_id, rownames(gh_meta))) |>
  select(-plant.id, -year, -swa, -pair, -row_id) |>
  as.matrix()
rownames(gh_raw_mat) <- rownames(gh_meta)

field_raw_mat <- field_wide |>
  mutate(row_id = paste(plant.id, year, collection.type, pair2, sep = "~")) |>
  filter(row_id %in% rownames(field_meta)) |>
  arrange(match(row_id, rownames(field_meta))) |>
  select(all_of(field_compounds)) |>
  as.matrix()
rownames(field_raw_mat) <- rownames(field_meta)

set.seed(4721)
gh_simper_raw    <- simper(gh_raw_mat,    group = gh_meta$swa,    permutations = 999)
field_simper_raw <- simper(field_raw_mat, group = field_meta$swa, permutations = 999)

cat("\n--- GH SIMPER (raw ng/g/hr, Bray-Curtis) ---\n")
print(summary(gh_simper_raw))

cat("\n--- Field SIMPER (raw ng/g/hr, Bray-Curtis) ---\n")
print(summary(field_simper_raw))

# --- 4B. CLR-based compound ranking (Aitchison analog) ---
# Ranks compounds by squared difference in mean CLR values between groups.
# Each compound's contribution to squared Aitchison distance = (mean_E- - mean_E+)^2,
# expressed as a percentage of the total squared centroid distance.
clr_rank <- function(clr_mat, meta_swa, label_absent, label_present) {
  mat          <- as.matrix(clr_mat)
  mean_absent  <- colMeans(mat[meta_swa == label_absent,  , drop = FALSE])
  mean_present <- colMeans(mat[meta_swa == label_present, , drop = FALSE])
  sq_diff      <- (mean_absent - mean_present)^2
  tibble(
    compound     = names(sq_diff),
    mean_absent  = round(mean_absent,  3),
    mean_present = round(mean_present, 3),
    sq_diff      = round(sq_diff,      4),
    pct_contrib  = round(sq_diff / sum(sq_diff) * 100, 2)
  ) |> arrange(desc(pct_contrib))
}

gh_clr_rank    <- clr_rank(gh_clr,    gh_meta$swa,    "Endophyte Absent", "Endophyte Present")
field_clr_rank <- clr_rank(field_clr, field_meta$swa, 0, 1)

cat("\n--- GH CLR compound ranking (Aitchison analog) ---\n")
print(gh_clr_rank, n = Inf)

cat("\n--- Field CLR compound ranking (Aitchison analog) ---\n")
print(field_clr_rank, n = Inf)


# ============================================================
# SECTION 5: INDICATOR VALUE (IndVal)
# ============================================================
# Identifies specific compounds significantly associated with
# E+ or E- groups. Useful for pinpointing candidate compounds
# for follow-up work.

gh_indval <- multipatt(
  as.data.frame(as.matrix(gh_clr)),
  gh_meta$swa,
  func    = "IndVal.g",
  control = how(nperm = 999)
)
cat("\n--- GH IndVal ---\n")
summary(gh_indval)

field_indval <- multipatt(
  as.data.frame(as.matrix(field_clr)),
  field_meta$swa,
  func    = "IndVal.g",
  control = how(nperm = 999)
)
cat("\n--- Field IndVal ---\n")
summary(field_indval)


# ============================================================
# SECTION 6: RESULTS SUMMARY
# ============================================================
# Helper to pull key stats from adonis2 output
extract_permanova_row <- function(aov_obj, row_term) {
  as.data.frame(aov_obj) |>
    rownames_to_column("term") |>
    filter(term == row_term) |>
    transmute(
      term,
      Df,
      F      = round(F, 3),
      R2_pct = round(R2 * 100, 2),
      p      = `Pr(>F)`
    )
}

gh_summary <- bind_rows(
  extract_permanova_row(gh_permanova, "swa"),
  extract_permanova_row(gh_permanova, "pair")
)

field_summary <- bind_rows(
  extract_permanova_row(field_permanova, "swa"),
  extract_permanova_row(field_permanova, "collection.type"),
  extract_permanova_row(field_permanova, "pair2")
)

cat("\n========================================\n")
cat("GREENHOUSE PERMANOVA (marginal, strata = year)\n")
cat("========================================\n")
print(gh_summary, row.names = FALSE)

cat("\n========================================\n")
cat("FIELD PERMANOVA (marginal, strata = year)\n")
cat("========================================\n")
print(field_summary, row.names = FALSE)

cat("\n--- GH betadisper (year / swa) ---\n")
cat("  year p =", betadisp_gh_year_p, "\n")
cat("  swa  p =", betadisp_gh_swa_p,  "\n")

cat("\n--- Field betadisper (year / type / swa) ---\n")
cat("  year p =", betadisp_field_year_p, "\n")
cat("  type p =", betadisp_field_type_p, "\n")
cat("  swa  p =", betadisp_field_swa_p,  "\n")
