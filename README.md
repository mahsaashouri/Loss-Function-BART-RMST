# Generalized Bayesian Additive Regression Trees for Restricted Mean Survival Time Inference

Prediction methods for time-to-event outcomes often utilize survival models that rely on strong assumptions about noninformative censoring or on how individual-level covariates and survival functions are related. When the main interest is in predicting individual-level restricted mean survival times (RMST), reliance on such assumptions can lead to poor predictive performance if these assumptions are not satisfied. We propose a generalized Bayes framework that avoids full probability modeling of all survival outcomes by using an RMST-targeted loss function that depends on a collection of inverse probability of censoring weights (IPCW). In our generalized Bayes formulation, we utilize a flexible additive tree regression model for the RMST function, and the posterior distribution of interest is obtained through model-averaging IPCW-conditional loss function-based pseudo-Bayesian posteriors. Because informative censoring can be captured by the IPCW-dependent loss function, our approach only requires one to specify a model for the censoring distribution, thereby obviating the need for complex joint modeling to handle informative censoring. We evaluate the performance of our method through a series of simulations that compare it with several well-known survival machine learning methods, and we illustrate the application of our method using a multi-site cohort of breast cancer patients with clinical and genomic covariates.

* Preprint is available [here](https://arxiv.org/abs/2402.17920)

## Dataset

Breast cancer outcomes assembled by the Molecular Taxonomy of Breast Cancer International Consortium (METABRIC) - available [here - accessed Jan 05, 2024](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric?resource=download)

## Reproducing results

* `NewFriedmanSimulationStudy_Independent.R`: results for the Friedman function simulation study for the case of noninformative censoring (*Tables 1*, *2*, and *4*).
* `ExpDependent_FriedmanSim.R`: results for the Friedman function simulation study for the case of informative censoring (*Table 3*).
* `NewMultinormAR1SimulationStudy.R`: results for the absolute value linear model with correlated predictors simulation study for the case of noninformative censoring (*Table 5*).
* `NEWMETABRIC.R`: results for the real-world application - METABRIC dataset (*Section 5*).


## Software

* We have implemented an R package `rmstbart` which is available for download [here](https://github.com/nchenderson/rmstbartold).
