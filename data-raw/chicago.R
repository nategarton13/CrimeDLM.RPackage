library("dplyr")

chicago <- readr::read_csv("crime_counts.csv") %>%
  mutate(type = tolower(`Primary Type`)) %>%
  select(year, month, type, count)

devtools::use_data(chicago, overwrite=TRUE)
