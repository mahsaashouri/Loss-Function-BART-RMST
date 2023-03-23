
library(tidyverse)
METABRIC <- read.csv('METABRIC_RNA_Mutation.csv', header = TRUE)[,c(1:31)]

METABRIC <- METABRIC %>% 
  select(-c('patient_id', 'cancer_type', 'cancer_type_detailed'))
## missing values in all columns are around 600 rows
METABRIC <- na.omit(METABRIC)

DATA <- METABRIC[ , !(names(METABRIC) %in% c('overall_survival'))]
DATA <- model.matrix(overall_survival_months~.-1, data = DATA)
Y <- METABRIC$overall_survival_months
delta <- METABRIC$overall_survival
sgrid <- seq(0, 4000, by=1)

tau <- 500
gam_alph <- 20
sigma <- 1.0

bcart_mod <- RMST_BCART(Y, delta, DATA, ndraws=500, tau=500, sigma.mu=1.2)
