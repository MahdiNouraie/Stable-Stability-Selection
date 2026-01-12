# This script generates Figure 4b in the paper 
# The script applies stable stability selection to the
# Riboflavin data using LASSO
# To reproduce the figure reported in the paper, 
# you need to set the number of
# subsamples B = 500. For faster computation, 
# you may reduce the number of subsamples.

options(warn=-1) # Turn off warnings
if (!requireNamespace("hdi")) {install.packages("hdi")} # install package if not already installed
if (!requireNamespace("glmnet")) {install.packages("glmnet")} # install package if not already installed
if (!requireNamespace("latex2exp")) {install.packages("latex2exp")} # install package if not already installed
if (!requireNamespace("ggplot2")) {install.packages("ggplot2")} # install package if not already installed


library(hdi) # Load the hdi package for riboflavin data
library(glmnet) # Load the glmnet package for LASSO
library(latex2exp) # Load the latex2exp package for TeX formatting
library(ggplot2) # Load the ggplot2 package for plotting

#sessionInfo()
#R version 4.3.2 (2023-10-31)
#Platform: aarch64-apple-darwin20 (64-bit)
#Running under: macOS 15.0.1

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: Australia/Sydney
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods  
#[7] base     

#other attached packages:
#  [1] ggplot2_3.5.1   latex2exp_0.9.6 glmnet_4.1-8    Matrix_1.6-1.1 
#[5] hdi_0.1-9       scalreg_1.0.1   lars_1.3       

#loaded via a namespace (and not attached):
#  [1] gtable_0.3.5      dplyr_1.1.4       compiler_4.3.2   
#[4] tidyselect_1.2.1  Rcpp_1.0.13       stringr_1.5.1    
#[7] parallel_4.3.2    splines_4.3.2     scales_1.3.0     
#[10] lattice_0.22-6    R6_2.5.1          linprog_0.9-4    
#[13] generics_0.1.3    shape_1.4.6.1     iterators_1.0.14 
#[16] MASS_7.3-60       tibble_3.2.1      munsell_0.5.1    
#[19] pillar_1.9.0      rlang_1.1.4       utf8_1.2.4       
#[22] stringi_1.8.4     cli_3.6.3         withr_3.0.1      
#[25] magrittr_2.0.3    foreach_1.5.2     grid_4.3.2       
#[28] rstudioapi_0.16.0 lifecycle_1.0.4   vctrs_0.6.5      
#[31] lpSolve_5.6.21    glue_1.8.0        codetools_0.2-20 
#[34] survival_3.7-0    fansi_1.0.6       colorspace_2.1-1 
#[37] tools_4.3.2       pkgconfig_2.0.3  

set.seed(26) # Set seed for reproducibility

data(riboflavin) # Load the riboflavin data
data <- as.data.frame(cbind(Y=riboflavin$y - 1, X=riboflavin$x)) # Convert the data to a data frame
rm(riboflavin) # Remove the original data to save memory
x <- as.matrix(data[, -1]) # Extract the predictors
y <- data[,1] # Extract the response

p <- ncol(x) # Number of predictors
B <- 500 # Number of subsamples
x <- scale(x) # Standardise the predictors
cv_lasso <- cv.glmnet(x, y, nfolds = 10, alpha = 1) # Perform cross-validation for LASSO

candidate_set <- cv_lasso$lambda # Candidate set of lambda values

S_list <- vector("list", length(candidate_set)) # Initialize a list to store the selection matrices
names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries

# Stability Selection for each lambda in candidate_set
for (lambda_idx in seq_along(candidate_set)) {
  
  lambda <- candidate_set[lambda_idx]  # Current lambda value
  S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
  colnames(S) <- colnames(x) # Set the column names of S to the predictors
  
  for (i in 1:B) {
    # Sub-sample the data (half of the original data without replacement)
    model_data <- data[sample(1:nrow(data), nrow(data) / 2, replace = FALSE), ]
    
    # Prepare the response and predictors
    x_sub <- as.matrix(model_data[, -1])  # Extract the predictors
    x_sub <- scale(x_sub) # Standardise the predictors
    y_sub <- model_data$Y  # Extract the response
    
    # Fit the LASSO model with the current lambda
    lasso_model <- glmnet(x_sub, y_sub, alpha = 1, lambda = lambda)
    
    # Extract significant predictors (ignoring the intercept, hence [-1])
    significant_predictors <- ifelse(coef(lasso_model) != 0, 1, 0)[-1]
    
    # Store the significant predictors in matrix S
    S[i, ] <- significant_predictors
  }
  
  # Store the matrix S for the current lambda in the corresponding list entry
  S_list[[lambda_idx]] <- S
}

# Stability measure (2018) from "https://github.com/nogueirs/JMLR2018/blob/master/R/getStability.R"
getStability <- function(X,alpha=0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability
  
  M<-nrow(X)
  d<-ncol(X)
  hatPF<-colMeans(X)
  kbar<-sum(hatPF)
  v_rand=(kbar/d)*(1-kbar/d)
  stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate
  
  ## then we compute the variance of the estimate
  ki<-rowSums(X)
  phi_i<-rep(0,M)
  for(i in 1:M){ 
    phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
  }
  phi_bar=mean(phi_i)
  var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate
  
  ## then we calculate lower and upper limits of the confidence intervals
  z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval
  
  return(list("stability"=stability,"variance"=var_stab,"lower"=lower,"upper"=upper))
  
}

stability_results <- lapply(S_list, function(S) { # Loop through regularisation set and compute stability values
  getStability(S) # Compute stability values
})


stab_values <- c() # Initialize an empty vector to store stability values of each lambda
for (i in 1:length(candidate_set)) {
  temp <- stability_results[[paste0("lambda_", i)]]$stability # Extract stability values
  stab_values <- c(stab_values, temp) # Append stability values to the vector
}

max_stability <- max(stab_values) # Find the maximum stability value
stability_1sd_threshold <- max_stability - sd(stab_values) # Define the stability threshold as max stability - 1SD
index_of_stable_1sd <- max(which(stab_values >= stability_1sd_threshold)) # since candidate_set is in decreasing order, 
#we find the index of the stable.1sd lambda value by maximum index

S_stable_1sd <- S_list[[index_of_stable_1sd]] # Extract the selection matrix for the stable.1sd lambda value
stability <- data.frame() # Initialize an empty data frame to store stability values
for (k in 2:nrow(S_stable_1sd)){ # Loop through sub-samples results for lambda stable.1sd
  output <- getStability(S_stable_1sd[1:k,]) # Compute stability values
  stability <- rbind(stability, data.frame(k, output$stability, output$variance, output$lower, output$upper)) # Append stability values to the data frame
}
colnames(stability) <- c('Iteration', 'Stability', 'Variance', 'Lower', 'Upper') # Set column names of the data frame

ggplot(stability, aes(x = Iteration, y = Stability)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = 'blue', alpha = 0.7) + # Add ribbon for confidence interval
  labs(title = TeX('Stability of Stability Selection ($\\lambda = \\lambda_{stable.1sd}$)'), 
       x = 'Iteration (sub-sample)', y = TeX('Stability ($\\hat{\\Phi}$)'))+
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),       # Title text size
    axis.title.x = element_text(size = 18),     # X-axis label size
    axis.title.y = element_text(size = 18),     # Y-axis label size
    axis.text.x = element_text(size = 16),      # X-axis tick text size
    axis.text.y = element_text(size = 16)       # Y-axis tick text size
  )
