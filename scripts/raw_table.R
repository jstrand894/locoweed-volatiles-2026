source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "00_load_data.R"))

joined.data

summary_table <- joined.data %>%
  filter(!is.na(swa)) %>%
  group_by(compound, swa) %>%
  summarise(
    mean = mean(ngghr, na.rm = TRUE),
    se = sd(ngghr, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(mean_se = sprintf("%.2f ± %.2f", mean, se)) %>%
  select(compound, swa, mean_se) %>%
  pivot_wider(
    names_from = swa,
    values_from = mean_se,
    names_prefix = "swa_"
  )

print(summary_table, n = 100)



# Prepare data for plotting
plot_data <- joined.data %>%
  filter(!is.na(swa)) %>%
  group_by(compound, swa) %>%
  summarise(
    mean = mean(ngghr, na.rm = TRUE),
    se = sd(ngghr, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(swa = factor(swa, labels = c("swa_0", "swa_1")))

# Create faceted bar chart
ggplot(plot_data, aes(x = swa, y = mean, fill = swa)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, linewidth = 0.8) +
  facet_wrap(~compound, scales = "free_y") +
  scale_fill_manual(values = c("swa_0" = "#2E86AB", "swa_1" = "#A23B72")) +
  labs(
    title = "ngghr by Compound and SWA Treatment",
    x = "SWA Treatment",
    y = "Mean ngghr",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )



field_table <- 
field_comb %>% 
  filter(!is.na(swa)) %>%
  group_by(compound, swa, collection.type) %>%
  summarise(
    mean = mean(ngghr, na.rm = TRUE),
    se = sd(ngghr, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(mean_se = sprintf("%.2f ± %.2f", mean, se)) %>%
  select(compound, swa, collection.type, mean_se) %>%
  pivot_wider(
    names_from = c(collection.type, swa),
    values_from = mean_se,
    names_sep = "_"
  )


field_comb %>%
  filter(!is.na(swa)) %>%
  filter(!is.na(compound)) %>%
  filter(compound == "z3hex") %>% 
  group_by(compound, swa, collection.type) %>%
  summarise(
    mean = mean(ngghr, na.rm = TRUE),
    se = sd(ngghr, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = compound, y = mean, fill = factor(swa))) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6,
           position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, linewidth = 0.8,
                position = position_dodge(width = 0.6)) +
  facet_wrap(~collection.type)
  
  
  
  ggplot(aes(x = compound, y = ngghr, fill = factor(swa))) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~collection.type) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Compound",
    y = "ngghr",
    fill = "SWA",
    title = "Volatile Compounds by Collection Type and SWA"
  )

write.csv(field_table, "field_table.csv", row.names = FALSE)






# looking into z3hex
# Identify the outlier


field_comb %>%
  filter(!(plant.id == 742 & compound == "z3hex" & ngghr > 300)) %>%
  filter(!is.na(swa)) %>%
  filter(!is.na(compound)) %>%
  filter(compound == "z3hex") %>% 
  group_by(compound, swa, collection.type) %>%
  summarise(
    n = n(),
    mean = mean(ngghr, na.rm = TRUE),
    se = sd(ngghr, na.rm = TRUE) / sqrt(n()),
    median = median(ngghr, na.rm = TRUE),
    sd = sd(ngghr, na.rm = TRUE),
    min = min(ngghr, na.rm = TRUE),
    max = max(ngghr, na.rm = TRUE),
    .groups = "drop"
  )

