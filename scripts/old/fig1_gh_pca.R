# ============================================================
# fig1_gh_pca.R
# Figure 1: Greenhouse PCA biplot (Spring 2022 + 2023)
#
# INPUTS (from 03_pca.R):
#   gh_pca_data, gh_pca_loadings, gh_pc1_var, gh_pc2_var
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "03_pca.R"))


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
    panel.grid        = element_blank(),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  )

GH_PCA_Plot

ggsave("figures/Fig1_GH_PCA.png", GH_PCA_Plot,
       width = 6, height = 5, dpi = 300)









GH_PCA_Plot <-
  ggplot(gh_pca_data) +
  geom_point(aes(PC1, PC2, shape = year, fill = swa),
             size = 2.5, stroke = 0.4, alpha = 0.8, color = "black") +
  stat_ellipse(aes(PC1, PC2, color = swa, group = interaction(swa, year)),
               level = 0.95) +
  scale_fill_manual(values  = c("Absent" = "#D55E00", "Present" = "#009E73")) +
  scale_color_manual(values = c("Absent" = "#D55E00", "Present" = "#009E73")) +
  scale_shape_manual(values = c("2022" = 21, "2023" = 24)) +
  scale_linetype_manual(values = c("2022" = "solid", "2023" = "dashed")) +
  labs(
    x        = paste0("PC1 (", gh_pc1_var, "%)"),
    y        = paste0("PC2 (", gh_pc2_var, "%)"),
    fill     = "Endophyte",
    color    = "Endophyte",
    linetype = "Year"
  ) +
  guides(color = guide_legend(override.aes = list(fill = c("#C1392B", "#2E86C1")))) +
  theme_bw() +
  theme(
    panel.grid        = element_blank(),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  )
GH_PCA_Plot
