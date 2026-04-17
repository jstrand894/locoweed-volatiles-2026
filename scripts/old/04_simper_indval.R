# ============================================================
# 04_simper_indval.R
# SIMPER, CLR-based compound ranking, and IndVal analysis.
#
# INPUTS (from 01_clr_processing.R):
#   gh_clr, gh_meta, gh_raw_mat
#   field_clr, field_meta, field_raw_mat
#
# OUTPUTS:
#   gh_simper_raw, field_simper_raw
#   gh_clr_rank, field_clr_rank
#   gh_indval, field_indval
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "01_clr_processing.R"))


# ============================================================
# SIMPER (raw ng/g/hr, Bray-Curtis)
# ============================================================
# CLR-transformed values can be negative, which breaks Bray-Curtis.
# Raw ng/g/hr values are non-negative and appropriate for simper().

set.seed(4721)
gh_simper_raw    <- simper(gh_raw_mat,    group = gh_meta$swa,    permutations = 999)
field_simper_raw <- simper(field_raw_mat, group = field_meta$swa, permutations = 999)

cat("\n--- GH SIMPER (raw ng/g/hr, Bray-Curtis) ---\n")
print(summary(gh_simper_raw))

cat("\n--- Field SIMPER (raw ng/g/hr, Bray-Curtis) ---\n")
print(summary(field_simper_raw))


# ============================================================
# CLR-based compound ranking (Aitchison analog of SIMPER)
# ============================================================
# Ranks compounds by squared difference in mean CLR values between groups.
# Each compound's contribution = (mean_E- − mean_E+)^2 as % of total
# squared centroid distance.

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
# IndVal — Indicator Value Analysis
# ============================================================
# Identifies compounds significantly associated with E+ or E- groups.
# Useful for pinpointing candidate compounds for follow-up.

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
