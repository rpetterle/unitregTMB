
## Read the original CSV file
periodontal <- read.csv("data-raw/DataSet_PD.csv")

periodontal$gender <- factor(periodontal$gender, 
                             levels = c(0, 1), 
                             labels = c("Male", "Female"))

periodontal$hba1c <- factor(periodontal$hba1c, 
                            levels = c(0, 1), 
                            labels = c("Controlled", "Uncontrolled"))

periodontal$smoker <- factor(periodontal$smoker, 
                             levels = c(0, 1), 
                             labels = c("Non-smoker", "Smoker"))

periodontal$tooth <- factor(periodontal$tooth, 
                            levels = c(1, 2, 3, 4), 
                            labels = c("Molar", "Premolar", "Canine", "Incisor"))

periodontal$tooth <- stats::relevel(periodontal$tooth, ref = "Canine")

usethis::use_data(periodontal, overwrite = TRUE)