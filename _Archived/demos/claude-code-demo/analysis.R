# Differential Expression Analysis
# R/4.4.0 on OSC Ascend

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# Load data
data <- read.csv("expression_data.csv")
cat("Data loaded:", nrow(data), "observations\n")
print(head(data))

# Descriptive statistics by group and gene
stats_summary <- data %>%
  group_by(gene, group) %>%
  summarise(
    n = n(),
    mean = mean(expression),
    sd = sd(expression),
    se = sd / sqrt(n),
    .groups = "drop"
  )

cat("\nDescriptive Statistics:\n")
print(stats_summary)

# Save summary statistics
write.csv(stats_summary, "results/summary_statistics.csv", row.names = FALSE)

# Differential expression t-tests for each gene
genes <- unique(data$gene)
de_results <- data.frame()

for (g in genes) {
  gene_data <- data %>% filter(gene == g)
  control <- gene_data %>% filter(group == "Control") %>% pull(expression)
  disease <- gene_data %>% filter(group == "Disease") %>% pull(expression)

  test <- t.test(disease, control)

  de_results <- rbind(de_results, data.frame(
    gene = g,
    control_mean = mean(control),
    disease_mean = mean(disease),
    log2FC = log2(mean(disease) / mean(control)),
    t_statistic = test$statistic,
    p_value = test$p.value,
    significant = test$p.value < 0.05
  ))
}

cat("\nDifferential Expression Results:\n")
print(de_results)

# Save DE results
write.csv(de_results, "results/differential_expression.csv", row.names = FALSE)

# Visualization 1: Box plots
p1 <- ggplot(data, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  facet_wrap(~gene, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Disease" = "#E64B35")) +
  theme_bw() +
  labs(
    title = "Gene Expression: Control vs Disease",
    x = "Group",
    y = "Expression Level"
  ) +
  theme(legend.position = "none")

ggsave("figures/boxplots.png", p1, width = 8, height = 4, dpi = 150)
cat("\nSaved: figures/boxplots.png\n")

# Visualization 2: Heatmap
data_wide <- data %>%
  select(gene, sample, expression) %>%
  pivot_wider(names_from = sample, values_from = expression)

# Convert to matrix with gene names as row names
mat <- as.matrix(data_wide[, -1])
rownames(mat) <- data_wide$gene

# Simple heatmap using base R
png("figures/heatmap.png", width = 600, height = 400)
heatmap(mat,
        scale = "row",
        main = "Expression Heatmap",
        col = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()
cat("Saved: figures/heatmap.png\n")

# Visualization 3: Bar plot with error bars
p3 <- ggplot(stats_summary, aes(x = gene, y = mean, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.8), width = 0.25) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Disease" = "#E64B35")) +
  theme_bw() +
  labs(
    title = "Mean Expression by Gene and Group",
    x = "Gene",
    y = "Expression Level (mean +/- SE)"
  )

ggsave("figures/barplot.png", p3, width = 6, height = 4, dpi = 150)
cat("Saved: figures/barplot.png\n")

cat("\nAnalysis complete!\n")
cat("Results saved to: results/\n")
cat("Figures saved to: figures/\n")

# Session info
cat("\nSession Info:\n")
sessionInfo()
