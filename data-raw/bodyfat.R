
## Loading extra packages
library(usethis)
library(dplyr)
library(tidyr)

raw_df <- read.csv("data-raw/qbmult_dataset.csv")

bodyfat <- raw_df %>%
  rename(
    id = ID, age = AGE, bmi = BMI, gender = SEX, ipaq = IPAQ,
    arms = ARMS, legs = LEGS, trunk = BODY, android = ANDROID, gynoid = GYNECOID
  ) %>%
  mutate(
    across(c(arms, legs, trunk, android, gynoid), ~ .x / 100),
    
    gender = factor(gender, levels = c(1, 2), labels = c("Female", "Male")),
    ipaq = factor(ipaq, 
                  levels = c(0, 1, 2), # Altere estes números conforme o seu CSV original
                  labels = c("Sedentary", "Insufficiently active", "Active"))
  )

bodyfat_long <- bodyfat %>%
  pivot_longer(
    cols = c(arms, legs, trunk, android, gynoid),
    names_to = "regions",
    values_to = "y"
  ) %>%
  mutate(
    regions = tools::toTitleCase(regions),
    regions = factor(regions, levels = c("Arms", "Legs", "Trunk", "Android", "Gynoid"))
  ) %>%
  arrange(id, regions)

usethis::use_data(bodyfat, overwrite = TRUE, compress = "xz")
usethis::use_data(bodyfat_long, overwrite = TRUE, compress = "xz")
