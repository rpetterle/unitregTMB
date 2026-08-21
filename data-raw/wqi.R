
## Loading extra packages
library(usethis)
library(dplyr)

raw_df <- read.csv("data-raw/DataSet_WQI.csv", header = TRUE)

wqi <- raw_df %>%
  mutate(
    id = factor(id),
    quarter = factor(quarter),
    
    location = factor(location, levels = c("Upstream", "Reservoir", "Downstream"))
  ) %>%
  arrange(id, quarter, location)

usethis::use_data(wqi, overwrite = TRUE, compress = "xz")
