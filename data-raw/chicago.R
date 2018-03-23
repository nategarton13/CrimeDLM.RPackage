library("dplyr")

chicago <- readr::read_csv("crime_counts.csv") %>%
  mutate(type = tolower(`Primary Type`),
         month = factor(month,
                        levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                   "Jul","Aug","Sep","Oct","Nov","Dec"))) %>%
  select(year, month, type, count)

devtools::use_data(chicago, overwrite=TRUE)
