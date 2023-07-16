
library(survival)
library(BART)
library(glmnet)
library(mboost)
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
install_github("nchenderson/AFTrees")
library(AFTrees)

## reading data
METABRIC <- read.csv('METABRIC_RNA_Mutation.csv', header = TRUE)[,c(1:41)]

######################
## data pre-processing
######################

## drop unwanted columns
METABRIC <- METABRIC %>% 
  select(-c('patient_id', 'cancer_type', 'cancer_type_detailed', 'her2_status_measured_by_snp6', 'er_status', 
            'death_from_cancer', 'X3.gene_classifier_subtype', 'oncotree_code', 'cellularity'))
## check number of missing in each column of the dataset
colSums(is.na(METABRIC))

## combine 4ER+ and 4ER- as 4 in integrative_cluster column
METABRIC$integrative_cluster[METABRIC$integrative_cluster %in% c("4ER+","4ER-")] <- "4"

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

## drop rows empty in type_of_surgery
table(METABRIC$type_of_breast_surgery) ## 16 rows empty
METABRIC <- subset(METABRIC, !(type_of_breast_surgery == ''))

## drop rows empty in er_status_measured_by_ihc
table(METABRIC$er_status_measured_by_ihc) ## 25 rows empty
METABRIC <- subset(METABRIC, !(er_status_measured_by_ihc == ''))

## make cohort as a factor column
METABRIC$cohort <- as.factor(METABRIC$cohort)

## flip 1 and 0 valves in overall survival
table(METABRIC$overall_survival)
METABRIC$overall_survival <- ifelse(METABRIC$overall_survival == 1, 0, ifelse(METABRIC$overall_survival == 0, 1, 
                                                                              METABRIC$overall_survival))
############################
## running different methods
############################

## split dataset into training and test set
set.seed(123)
## 400 rows as test set (~20%)
DATA <- METABRIC[ , !(names(METABRIC) %in% c('overall_survival'))]
DATA <- model.matrix(overall_survival_months~.-1, data = DATA)


# Define the number of iterations and the proportion of data to be used for training
n_iterations <- 1
train_prop <- 0.7
sgrid <- seq(0, 4000, by=1)
tau <- 50
gam_alph <- 20
sigma <- 1.0
ndraws <- 500

bcart_fitted <- bart_fitted <- AFT_BART_fitted <- list()
# Loop through the iterations
for (i in 1:n_iterations) {
  
  # Randomly select index numbers for the training and test sets
  train_idx <- sample(1:nrow(DATA), round(train_prop * nrow(DATA)), replace = FALSE)
  test_idx <- setdiff(1:nrow(DATA), train_idx)
  
  # Create training and test sets using the selected index numbers
  train.set <- DATA[train_idx, ]
  test.set <- DATA[test_idx, ]
  
  Y <- METABRIC[train_idx,]$overall_survival_months
  delta <- METABRIC[train_idx,]$overall_survival
  
  Y.test <- METABRIC[test_idx,]$overall_survival_months
  #mu.test <- log(Y.test)
  mu.test <- (Y.test/gam_alph)*pgamma(tau, shape = gam_alph+1, rate = Y.test) + 
    tau*pgamma(tau, shape = gam_alph, rate = Y.test, lower.tail = FALSE)
  delta.test <- METABRIC[test_idx,]$overall_survival
  
  ## BCART
  #bcart_mod <- RMST_BCART(Y, delta, train.set, test.set,ndraws=ndraws, tau=tau)
  ## If doing dependent censoring use:
  bcart_mod <- RMST_BCART(Y, delta, train.set, test.set,
                          ndraws=ndraws, ipcw="dependent", tau=tau)
  #bcart_fitted[[i]] <- rowMeans(bcart_mod$fitted.values.test)
  bcart_fitted[[i]] <- bcart_mod
  
  ## BART
  bart_mod_ind <- RMST_BART(Y, delta, train.set, test.set, ndraws=ndraws, tau=tau, ntrees = 200)
  ## If doing dependent censoring use:
  bart_mod_dep <- RMST_BART(Y, delta, train.set, test.set, ntrees = 200,
                        ndraws=ndraws, ipcw="dependent", tau=tau)
  #bart_fitted[[i]] <- rowMeans(bart_mod$fitted.values.test)
  bart_fitted[[i]] <- bart_mod
  
  #### 3. AFT_BART model
  AFT_BART <- abart(train.set, Y, delta, x.test=test.set)
  ndraw_abart <- nrow(AFT_BART$yhat.test)
  AFT_fit_reps <- matrix(NA, nrow=ndraw_abart, ncol=nrow(test.set))
  for(k in 1:ndraw_abart) {
    aft_bart_mu <- AFT_BART$yhat.test[k,]
    aft_bart_sig <- AFT_BART$sigma[k]
    aft_bart_sigsq <- aft_bart_sig*aft_bart_sig
    gt_prob <- pnorm((log(tau) - aft_bart_mu)/aft_bart_sig, lower.tail=FALSE)
    lt_prob <- pnorm((log(tau) - aft_bart_sigsq - aft_bart_mu)/aft_bart_sig)
    
    AFT_fit_reps[k,] <- exp(aft_bart_sigsq/2 + aft_bart_mu)*lt_prob + tau*gt_prob
  }
  #AFT_BART_fitted[[i]] <- colMeans(AFT_fit_reps)
  AFT_BART_fitted[[i]] <- AFT_fit_reps
}


## plotting the first 10 repeated variables in one iteration - BART
VarImp <- tail(sort(colSums(bart_fitted[[1]]$split.vars)),10)

library(ggplot2)
# Create a data frame with the numbers and names
VarImpDataF <- data.frame(numbers = c(VarImp)/ndraws,
                 names = c(names(VarImp)))
library(ggrepel)
# Create the plot
ggplot(VarImpDataF, aes(x = seq_along(numbers), y = numbers)) +
  geom_line(size = 1, color = "black") +
  geom_point(size = 3, color = "blue") +
  #geom_text(aes(label = names), hjust = -.5, vjust = -.9, size = 6, color = "blue") +
  #geom_text_repel(aes(label = names), size = 6, color = "blue") +
  geom_label_repel(aes(label = names), size = 6, color = "blue") +
  labs(x = "Variable", y = "Count") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title=element_text(size=15), 
    legend.text=element_text(size=15))

## Average RMST
bcart.r <- rowMeans(bcart_mod$fitted.values.test)
sqrt(mean((bcart.r - mu.test)*(bcart.r - mu.test)))
bart.r <- rowMeans(bart_mod$fitted.values.test)
sqrt(mean((bart.r - mu.test)*(bart.r - mu.test)))
AFT.r <- colMeans(AFT_fit_reps)
sqrt(mean((AFT.r - mu.test)*(AFT.r - mu.test)))

## Error 
delta_alpha <- 1
Gvec.test <- DrawIPCW(U=Y.test, delta=delta.test, Utau=pmin(Y.test, tau), sgrid=sgrid,
                 kappa0=1, delta_alpha=delta_alpha)
## BART
mean.test.bart <- rowMeans(bart_fitted[[1]]$fitted.values.test)

# calculate IPCW test performance - BART # 1151.15
IPCW.test.bart <- sum((delta.test/Gvec.test)*(Y.test-mean.test.bart)^2)/length(Y.test)

## BACRT
mean.test.bcart <- rowMeans(bcart_fitted[[1]]$fitted.values.test)

# calculate IPCW test performance - BCART # 948.4491
IPCW.test.bcart <- sum((delta.test/Gvec.test)*(Y.test-mean.test.bcart)^2)/length(Y.test)

## Confidence interval for each case - plot 
means <- rowMeans(bart_fitted[[1]]$fitted.values.test)
ses <- matrixStats::rowSds(bart_fitted[[1]]$fitted.values.test, na.rm=TRUE)
df_summary <- data.frame(row = 1:nrow(bart_fitted[[1]]$fitted.values.test), mean = means, se = ses)

df_summary$max <- df_summary$mean + 1.96 * df_summary$se
df_summary$min <- df_summary$mean - 1.96 * df_summary$se
df_summary <- subset(df_summary, select = -c(se))
## plot sample of lines
sample_idx <- sample(1:nrow(df_summary), round(0.01 * nrow(df_summary)), replace = FALSE)
sample_df_summary <- df_summary[sample_idx,]
#sample_long <- reshape2::melt(sample_df_summary, id.vars = "row")

# Sort the rows of df_summary by the mean column
df_summary_sorted <- df_summary[order(df_summary$mean),]

# create the plot
plot(1, type='n', xlim=c(min(df_summary_sorted$min), max(df_summary_sorted$max)), 
     ylim=c(0.5, nrow(df_summary_sorted)+0.5), ylab='Case number', xlab='Test set credible intervals for RMST')

# loop through the rows and draw a line for each
for (i in 1:nrow(df_summary)) {
  segments(df_summary_sorted$min[i], i, df_summary_sorted$max[i], i)
  points(df_summary_sorted$mean[i], i, pch=19, col='blue')
}

## Partial Dependence plots

Col_ParDep <- c('nf1','nbn', 'age_at_diagnosis', 'nottingham_prognostic_index') 
DATA <- METABRIC[ , !(names(METABRIC) %in% c('overall_survival'))]
DATA <- model.matrix(overall_survival_months~.-1, data = DATA)
ff <- matrix(NA, nrow = 10, ncol = 3)
partial_results <- list()
for(i in 1:length(Col_ParDep)){
  for(k in 1:10){
    pp <- seq(min(METABRIC[,Col_ParDep[i]]),max(METABRIC[,Col_ParDep[i]]), length.out = 10)
    Xtmp <- DATA
    Xtmp[,Col_ParDep[i]] <- rep(pp[k], nrow(DATA))
    Y <- METABRIC$overall_survival_months
    delta <- METABRIC$overall_survival
    bart_mod <- RMST_BART(Y, delta, Xtmp, tau=tau, ntrees = 200)
    ff[k,1] <- mean(rowMeans(bart_mod$fitted.values))
    ff[k,2] <- pp
    ff[k,3] <- Col_ParDep[i]
  } 
  partial_results[[length(partial_results)+1]] <- ff 
}

partial_results <- do.call("rbind", partial_results)

library(ggplot2)
library(patchwork)

# Create a list to store the individual plots
plots <- list()

# Iterate over each category and create a plot
for (category in unique(partial_results$index)) {
  subset_data <- subset(partial_results, index == category)
  
  plot <- ggplot(subset_data, aes(x = Value, y = MeanPrediction)) +
    geom_line() +
    xlab(category)+
    theme(axis.title = element_text(size = 22),  # Adjust the size of the axis titles
          axis.text = element_text(size = 20)) 
  
  # Set y-axis label for the leftmost plots
  if (category %in% unique(partial_results$index)[c(1,3)]) {
    plot <- plot + ylab("Mean Prediction") #+ theme(axis.title.y = element_blank())
  } else {
    plot <- plot + theme(axis.title.y = element_blank())
  }
  
  plots[[category]] <- plot
}

# Combine the individual plots using patchwork
combined_plot <- wrap_plots(plots, ncol = 2)

# Display the combined plot
combined_plot


## Posterior mean comparison

PostMean_ind <- c(rowMeans(bart_mod_ind$fitted.values.test))
PostMean_dep <-  c(rowMeans(bart_mod_dep$fitted.values.test))
PostMean_data <- data.frame('Independent' = PostMean_ind, 'Dependent' = PostMean_dep)
#write.csv(PostMean_data, 'PostMean_data_METABRIC.train.csv')

# Plot the data

ggplot(PostMean_data, aes(x = Independent, y = Dependent)) +
  geom_point() +
  labs(x = "Independent", y = "Dependent", title = "") +
  theme(axis.title = element_text(size = 14),  # Adjust the size of the axis titles
        axis.text = element_text(size = 12))   # Adjust the size of the axis labels
## Extra codes

#ff_1 <- c(44.64780, 44.26682, 44.55446, 44.53387, 44.25460, 44.04177, 43.84110, 44.01009, 44.54170, 44.09558)
#pp_1 <- c(-2.96860000, -1.95224444 ,-0.93588889,  0.08046667 , 1.09682222 , 2.11317778 , 3.12953333  ,4.14588889 , 5.16224444,  6.17860000)
#i_1 <- rep('nf1', 10)
#partial_results_1 <- data.frame(Value = pp_1,
#                                MeanPrediction = ff_1, index = i_1)

#ff_2 <-  c(44.54398, 44.25755, 44.38472, 44.22499, 44.41672, 44.20168, 44.33344, 44.20899, 44.14735, 44.48058)
#pp_2 <- c(-3.6898000, -2.6830111, -1.6762222, -0.6694333,  0.3373556,  1.3441444,  2.3509333,  3.3577222,  4.3645111,  5.3713000)
#i_2 <- rep('nbn', 10)
#partial_results_2 <- data.frame(Value = pp_2,
#                                MeanPrediction = ff_2, index = i_2)

#ff_3 <- c(44.02239, 44.25650, 44.65303, 44.31947, 44.27197, 44.14585, 44.17240, 44.08893, 44.01477, 44.04726)
#pp_3 <-  c(21.93000, 30.19222, 38.45444, 46.71667, 54.97889, 63.24111, 71.50333, 79.76556, 88.02778, 96.29000)
#i_3 <- rep('age_at_diagnosis', 10)
#partial_results_3 <- data.frame(Value = pp_3,
#                                MeanPrediction = ff_3, index = i_3)

#ff_4 <-  c(44.34663, 44.28516, 44.07619, 44.16539, 44.31767, 44.24659, 44.46195, 44.13971, 44.62427, 44.35170)
#pp_4 <-   c(1.020000, 1.613333, 2.206667, 2.800000, 3.393333, 3.986667, 4.580000, 5.173333, 5.766667, 6.360000)
#i_4 <- rep('nottingham_prognostic_index', 10)
#partial_results_4 <- data.frame(Value = pp_4,
#                               MeanPrediction = ff_4, index = i_4)

#partial_results <- rbind(partial_results_1, partial_results_2, partial_results_3, partial_results_4)

