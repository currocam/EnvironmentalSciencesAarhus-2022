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
  return(tibble::as_tibble(df))
}
