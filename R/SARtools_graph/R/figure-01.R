# Load libraries
library(tidyverse)
# Load custom function read_SARtools_into_tibble
source("SARtools_graph/R/common.R")
source("SARtools_graph/R/figure-01_colors.R")

# Configuration
# Save dir path into variable
## Edit this variable pls
dirpath <- "SARtools_graph/data/"

# Define custom aesthetics
axis_text_size_in_pt <- 15 
axis_title_size_in_pt <- 20
legend_text_size_in_pt <- 10 
n_rows_in_wrap <- 6
final_plot_scale <- 5
final_plot_width_in_inches <- 4
final_plot_height_in_inches <- 2
final_plot_dpi <- "retina"
final_plot_name_without_extension <- "SARtools_graph/figures/alfalfa_plot"

# Reading data
## Get files
txt_files <- list.files(dirpath,pattern = "*.txt")
## Read all files into long table
data <- paste0(dirpath, txt_files) %>%
  magrittr::set_names(
    txt_files %>%
      stringr::str_sub(1, 2)
    ) %>%
  purrr::map_dfr(
    .id = "file",
    read_SARtools_into_tibble)

# Check before plot
if(length(levels_colors) != length(levels(data$Lvl_2_letter))){
  stop(
    "There are a number of different levels and colors. Consider adding ",
    "or removing colors or coloring automatically by removing calls to the ",
    "scale_color_manual function.")   
}
# Plotting figures
## Save ggplot into figure variable
figure <- data %>%
  group_by(COG_Category_lvl_1, COG_Category_lvl_2_2,
           Up_or_down_regulated, file) %>%
  summarise(
    Frequency = sum(Count),
    Lvl_2_letter = unique(Lvl_2_letter)
  ) %>%
  mutate(
    Frequency = ifelse(Up_or_down_regulated == "Up", Frequency, -Frequency)
  ) %>%
  # Filtering [S] function out
  filter(Lvl_2_letter != "[S]") %>%
  ggplot(aes(x = Lvl_2_letter, y = Frequency)) +
  geom_bar(aes(fill = COG_Category_lvl_2_2, color = COG_Category_lvl_2_2),
           stat = "identity", position = "identity"
  ) +
  ylab("Frequency of up/down regulated genes") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_text(aes(y = Frequency, label = abs(Frequency))) +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values  = levels_colors) +
  scale_color_manual(values  = levels_colors) +
  xlab("Category level 2") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = axis_text_size_in_pt),
    axis.title.y = element_text(size = axis_title_size_in_pt),
    axis.title.x = element_text(size = axis_title_size_in_pt),
    legend.text = element_text(size = legend_text_size_in_pt),
  )+
  facet_wrap(~file) +
  guides(fill = guide_legend(nrow = n_rows_in_wrap, byrow = TRUE))

## Save ggplot in several formats
c("png", "svg", "pdf") %>%
  walk(
    ~ ggsave(
      paste(final_plot_name_without_extension, .x,sep = "."),
      plot = figure,
      device = .x,
      dpi = final_plot_dpi,
      scale = final_plot_scale,
      width = final_plot_width_in_inches,
      height = final_plot_height_in_inches
    )
  )

  