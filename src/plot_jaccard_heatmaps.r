library(readr)
library(dplyr)
library(ggplot2)

data_dir  <- "output"

tissue_a  <- "Artery_aorta" #tissue 1
tissue_b  <- "Artery_coronary"               # set tissue 2
resolutions <- c(1.5)   # set the resolutions you want

out_path <- file.path(data_dir, paste0(tissue_a, "_", tissue_b, "_jaccard_heatmap_", resolutions[1], ".pdf"))

all_df <- list()

for (res in resolutions) {
  path <- file.path(data_dir, paste0("community_comparisons_condor_res", res, ".tsv"))
  if (!file.exists(path)) {
    message("missing file: ", path)
    next
  }
  df_res <- read_tsv(path, show_col_types = FALSE) %>%
    mutate(resolution = res)
  all_df[[length(all_df) + 1]] <- df_res
}

df <- bind_rows(all_df)

# subset to the two tissues
df_sub <- df %>%
  filter(tissue1 %in% c(tissue_a, tissue_b),
         tissue2 %in% c(tissue_a, tissue_b))

# community+resolution labels for each side
df_sub <- df_sub %>%
  mutate(
    label1 = paste(tissue1, community1, resolution, sep = ":"),
    label2 = paste(tissue2, community2, resolution, sep = ":"),
    # canonical pair id: order-independent
    pair_id = if_else(label1 <= label2,
                      paste(label1, label2, sep = " | "),
                      paste(label2, label1, sep = " | "))
  )

# no aggregation: each row is a community pairing at a given resolution
df_plot <- df_sub

# order rows by resolution, then by jaccard (high to low) within resolution
row_order <- df_plot %>%
  arrange(resolution, desc(jaccard_index)) %>%
  pull(pair_id) %>%
  unique()

df_plot <- df_plot %>%
  mutate(
    pair_id  = factor(pair_id, levels = row_order),
    dummy_x  = factor("one")     # single x position → effectively one-axis plot
  )

n_high <- sum(df_plot$jaccard_index > 0.5, na.rm = TRUE)

p <- ggplot(df_plot, aes(x = 1, y = pair_id, fill = jaccard_index)) +
  geom_tile(color = NA) +   # no grid/borders
  scale_fill_gradient(
    low    = "white",
    high   = "darkred",
    limits = c(0, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.title.y  = element_text(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.x  = element_blank(),
    panel.grid    = element_blank()
  ) +
  labs(
    title    = paste("Jaccard per community pairing for", tissue_a, "and", tissue_b,
                     "at resolution", resolutions[1]),
    subtitle = paste("Pairs with Jaccard >", 0.5, ":", n_high),
    y        = "Community pairing",
    fill     = "Jaccard"
  )
ggsave(out_path, p, width = 10, height = 10, dpi = 300)