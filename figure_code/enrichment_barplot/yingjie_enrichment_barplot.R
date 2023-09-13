library(tidyverse)
setwd('C:/Users/flyku/Desktop/yingjie_plot/enrichment_barplot')
set.seed(42)


# Process the data
df1 <- df1 %>%
  slice(1:10) %>%
  mutate(group = 'up')

df2 <- df2 %>%
  slice(1:10) %>%
  mutate(group = 'down')

new_df <- rbind(df1, df2) %>%
  mutate(Term = str_replace_all(Term, " \\(GO.*", "")) %>%
  mutate(Adjusted.P.value = -1*log10(Adjusted.P.value),
         group = as.factor(group),
         Term = fct_inorder(Term)) %>%
  arrange(group, desc(Term))  # Rearrange the rows based on the 'group' and 'Term' columns

# Reverse the order of the rows
new_df$Term <- factor(new_df$Term, levels = rev(levels(new_df$Term)))

# Create the barplot
p2 <- ggplot(new_df,
             aes(x = Adjusted.P.value, y = Term, fill = group)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("") +
  labs(size = "Overlapping count") +
  theme(
    legend.title = element_text(size = 18, ),
    legend.text = element_text(size = 14, ),
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 18)
  ) +
  scale_x_continuous(name = "Combined score") +
  scale_size(range = c(8, 14))

# Plot
p2


# save plot

png(
  paste0('./yingjie_enrichment_barplot.png'),
  width = 4000,
  height = 3500,
  res = 300
)

p2
dev.off()

