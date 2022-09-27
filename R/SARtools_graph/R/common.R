# Define custom function for reading SARtools file
read_SARtools_into_tibble <- function(path){
  # Columns to select
  colnames_of_interest <- c(
    "COG_Category_lvl_1", "COG_Category_lvl_2_2", "Lvl_2_letter",
    "COG_Category_lvl_3", "Count", "Up_or_down_regulated"
  )
  # Read tsv with proper col types
  df <- read_tsv(
    path,
    col_types = c(
      .default = "factor",
      Count = "numeric"
    )
  ) %>%
    select(colnames_of_interest)
  # Arrange factors
  levels(df$Up_or_down_regulated) <- c("Up", "Down") # Show first up, then Down
  df$COG_Category_lvl_1 <- fct_infreq(df$COG_Category_lvl_1) #Show most frequent first
  df$COG_Category_lvl_2_2 <- fct_reorder(
    fct_infreq(df$COG_Category_lvl_2_2), # Show most frequent first
    as.numeric(df$COG_Category_lvl_1) # taking in account category level 1
  )
  df$Lvl_2_letter <- fct_reorder( # Again the same
    fct_infreq(df$Lvl_2_letter), as.numeric(df$COG_Category_lvl_1)
  )
  return(tibble::as_tibble(df))
}
