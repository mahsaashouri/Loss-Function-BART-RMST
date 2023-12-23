# Generalized Bayesian Additive Regression Trees for Restricted Mean Survival Time Prediction

Prediction methods for time-to-event outcomes often utilize 
survival models that
rely on strong assumptions about noninformative censoring 
or on how individual-level covariates
and survival functions are related. 
When the main interest is in predicting an individual-level expected
time to failure, reliance on such assumptions can lead
to poor predictive performance if these assumptions
are not satisfied.
By using a loss function that depends on the inverse probability of censoring weights (IPCW) instead of using a full probability model for all survival outcomes, we propose a generalized Bayes framework that uses an additive regression tree model
that explicitly targets predicting the restricted mean survival time. 
In our generalized Bayes formulation, the posterior distribution of interest
is obtained through model-averaging of IPCW-conditional
loss function-based Gibbs posteriors where the model weights
are based on a Gamma process model for the censoring distributions. 
Because informative censoring is captured by the IPCW-dependent
loss function, our approach only requires one to specify a
model for the censoring distribution, thereby obviating 
the need for complex joint modeling to handle informative censoring.
We evaluate the performance of our method through a series of simulations 
and compare its performance with several well-known survival machine learning methods. We illustrate the application of our method using a
mult-site cohort of breast cancer patients that has both clinical 
and genomic covariates.

* We have implemented an R package which is available for download at \url{github.com/nchenderson/BARTTrial}.
