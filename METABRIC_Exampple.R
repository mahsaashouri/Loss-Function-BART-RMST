
library(tidyverse)
setwd("~/Documents/LossFunctionBART/Loss-Function--BART")
source("A_matrix.R")
source("ChangeMove.R")
source("D_probabilities.R")
source("GPDraw.R")
source("GrowMove.R")
source("LogLik.R")
source("PruneMove.R")
#source("PruneMove_V2.R")
source("RMST_BCART.R")
source("RMST_BART.R")
source("RMST_MHRatio.R")
source("Recursive_A_matrix.R")
source("Transition_Prob.R")
source("fitted_value.R")
source("prior_conditional_on_D_V2.R")
source("tree-configuration.R")
source("DrawIPCW.R")

## reading data
METABRIC <- read.csv('METABRIC_RNA_Mutation.csv', header = TRUE)[,c(1:31)]

######################
## data pre-processing
######################

## drop unwanted columns
METABRIC <- METABRIC %>% 
  select(-c('patient_id', 'cancer_type', 'cancer_type_detailed', 'her2_status_measured_by_snp6', 'er_status_measured_by_ihc', 
            'death_from_cancer', 'X3.gene_classifier_subtype'))
## check number of missing in each column of the dataset
colSums(is.na(METABRIC))

## drop tumor_stage = 0 
table(METABRIC$tumor_stage)
METABRIC <- subset(METABRIC, !(tumor_stage == 0 & !is.na(tumor_stage)))

## replace tumor_stage = NA with cannot be evaluated
METABRIC$tumor_stage[is.na(METABRIC$tumor_stage)] <- 0
## make tumor_stage as a factor column
METABRIC$tumor_stage <- as.factor(METABRIC$tumor_stage)

## drop rows with missing tumor size
table(METABRIC$tumor_size)
METABRIC <- subset(METABRIC, !is.na(tumor_size))

## replace neoplasm = NA with 0
table(METABRIC$neoplasm_histologic_grade)
METABRIC$neoplasm_histologic_grade[is.na(METABRIC$neoplasm_histologic_grade)] <- 0
## make neoplasm as a factor column
METABRIC$neoplasm_histologic_grade <- as.factor(METABRIC$neoplasm_histologic_grade)

## replace mutation_count = NA with 0
table(METABRIC$mutation_count)
METABRIC$mutation_count[is.na(METABRIC$mutation_count)] <- 0


## running the analysis

DATA <- METABRIC[ , !(names(METABRIC) %in% c('overall_survival'))]
DATA <- model.matrix(overall_survival_months~.-1, data = DATA)
Y <- METABRIC$overall_survival_months
delta <- METABRIC$overall_survival
sgrid <- seq(0, 4000, by=1)

tau <- 500
gam_alph <- 20
sigma <- 1.0

bcart_mod <- RMST_BCART(Y, delta, DATA, ndraws=500, tau=500, sigma.mu=1.2)
