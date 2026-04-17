# ============================================================
# fig2_field_pca.R
# Figure 2: Field PCA biplot (2022 + 2023), faceted by tissue type
#
# INPUTS (from 03_pca.R):
#   field_pca_data, field_pca_loadings, field_pc1_var, field_pc2_var
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "03_pca.R"))


FIELD_PCA_Plot <-
  ggplot(field_pca_data) +
  geom_point(aes(PC1, PC2, shape = year, fill = swa),
             size = 2.5, stroke = 0.4, alpha = 0.8, color = "black") +
  stat_ellipse(aes(PC1, PC2, color = swa, group = swa),
               level = 0.95) +
  scale_fill_manual(values  = c("Absent" = "#D55E00", "Present" = "#009E73")) +
  scale_color_manual(values = c("Absent" = "#D55E00", "Present" = "#009E73")) +
  scale_shape_manual(values = c("2022" = 21, "2023" = 24)) +
  facet_wrap(~ collection.type) +
  labs(
    x     = paste0("PC1 (", field_pc1_var, "%)"),
    y     = paste0("PC2 (", field_pc2_var, "%)"),
    fill  = "Endophyte",
    color = "Endophyte",
    shape = "Year"
  ) +
  guides(
    shape = guide_legend(override.aes = list(fill = "gray40")),
    color = guide_legend(override.aes = list(fill = c("#D55E00", "#009E73")))
  ) +
  theme_bw() +
  theme(
    panel.grid        = element_blank(),
    strip.background  = element_rect(fill = NA),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  )

FIELD_PCA_Plot

ggsave("figures/Fig2_Field_PCA.png", FIELD_PCA_Plot,
       width = 8, height = 5, dpi = 300)
