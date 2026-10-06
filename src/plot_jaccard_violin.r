rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(stringr)
library(tidyr)

data_dir <- "output"

# tissues you want to compare (edit as needed)
tissues <- c("Artery_aorta", "Artery_coronary", "Brain_cerebellum")
resolutions <- c(1, 1.25, 1.5, 1.75, 2, 5)

out_path <- file.path(
  data_dir,
  paste0(paste(tissues, collapse = "_"), "_jaccard_violin_all_pairs_all_resolutions.pdf")
)

# read all requested resolutions and tag with resolution value
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

# cross-tissue pairs among the chosen tissues (no within-tissue)
df_pairs <- df %>%
  filter(
    tissue1 %in% tissues,
    tissue2 %in% tissues,
    tissue1 != tissue2
  ) %>%
  mutate(
    pair = if_else(tissue1 <= tissue2,
                   paste(tissue1, tissue2, sep = " vs "),
                   paste(tissue2, tissue1, sep = " vs "))
  )

# p <- ggplot(df_pairs, aes(x = pair, y = jaccard_index)) +
#   geom_violin(fill = "grey80", color = "black", trim = FALSE) +
#   theme_minimal() +
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor   = element_blank()
#   ) +
#   scale_y_continuous(limits = c(0, 1)) +          # full range
#   scale_y_break(c(0.15, 1)) +                    # <-- interval to remove
#   labs(
#     title = paste("Jaccard distributions across resolutions for tissue pairs:",
#                   paste(tissues, collapse = ", ")),
#     x     = "Tissue pair",
#     y     = "Jaccard index"
#   )

# ggsave(out_path, p, width = 7, height = 4, dpi = 300)

p_density_pairs <- ggplot(df_pairs, aes(x = jaccard_index, color = pair)) +
  geom_density(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = paste(
      "Density of Jaccard indices per tissue pair across resolutions:",
      paste(tissues, collapse = ", ")
    ),
    x = "Jaccard index",
    y = "Density",
    color = "Tissue pair"
  )

ggsave(
  file.path(
    data_dir,
    paste0(paste(tissues, collapse = "_"), "_jaccard_density_by_pair_all_resolutions.pdf")
  ),
  p_density_pairs,
  width = 7,
  height = 4,
  dpi = 300
)

# keep only Jaccard indices >= 0.2
df_pairs_filt <- df_pairs %>%
  filter(jaccard_index >= 0.2)

# boxplot instead of violin
p_box <- ggplot(df_pairs_filt, aes(x = pair, y = jaccard_index)) +
  geom_boxplot(fill = "grey80", color = "black") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = paste(
      "Jaccard boxplots (J ≥ 0.2) across resolutions"),
    x = "Tissue pair",
    y = "Jaccard index"
  )

ggsave(
  file.path(
    data_dir,
    paste0(paste(tissues, collapse = "_"), "_jaccard_boxplot_all_pairs_all_resolutions_Jge0.2.pdf")
  ),
  p_box,
  width = 7,
  height = 4,
  dpi = 300
)