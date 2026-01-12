# This script generates Figure 4a in the paper.
# The script generates synthetic data with 500 predictors, 50 observations, and 
#a correlation coefficient of 0.2 followed by stable stability selection 
# To reproduce the figure reported in the paper, you need to set the number of
# subsamples B = 500. For faster computation, you may reduce the number of subsamples.

options(warn=-1) # Turn off warnings
if (!requireNamespace("MASS")) {install.packages("MASS")} # install package if not already installed
if (!requireNamespace("glmnet")) {install.packages("glmnet")} # install package if not already installed
if (!requireNamespace("latex2exp")) {install.packages("latex2exp")} # install package if not already installed
if (!requireNamespace("ggplot2")) {install.packages("ggplot2")} # install package if not already installed

# Load necessary library
library(MASS) # load the MASS package for mvrnorm
library(glmnet) # load the glmnet package for LASSO
library(latex2exp) # load the latex2exp package for TeX
library(ggplot2) # load the ggplot2 package for plotting

# sessionInfo() # to check the version of R and the packages
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
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#other attached packages:
#  [1] ggplot2_3.5.1   latex2exp_0.9.6 glmnet_4.1-8    Matrix_1.6-1.1  MASS_7.3-60    
#loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5       cli_3.6.3         rlang_1.1.4       stringi_1.8.4    
#[5] generics_0.1.3    glue_1.8.0        colorspace_2.1-1  fansi_1.0.6      
#[9] scales_1.3.0      grid_4.3.2        tibble_3.2.1      munsell_0.5.1    
#[13] foreach_1.5.2     lifecycle_1.0.4   stringr_1.5.1     compiler_4.3.2   
#[17] dplyr_1.1.4       codetools_0.2-20  pkgconfig_2.0.3   Rcpp_1.0.13      
#[21] rstudioapi_0.16.0 lattice_0.22-6    R6_2.5.1          tidyselect_1.2.1 
#[25] utf8_1.2.4        pillar_1.9.0      splines_4.3.2     shape_1.4.6.1    
#[29] magrittr_2.0.3    withr_3.0.1       tools_4.3.2       gtable_0.3.5     
#[33] iterators_1.0.14  survival_3.7-0  

set.seed(26) # Set seed for reproducibility
n <- 50 # Number of samples
p <- 500 # Number of predictors
rho <- 0.2 # correlation coefficient

# Function to generate synthetic data
generator <- function(n, p, rho){
  # Construct the covariance matrix with rho^(|i-j|) correlation
  Sigma <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- rho^abs(i - j)
    }
  }
  predictors <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma) # Generate multivariate normal predictors
  beta <- c(1.5, 1.1, rep(0, p - 2)) # True regression coefficients
  eps <- rnorm(n, mean = 0, sd = 1) # Error term
  response <- predictors %*% beta + eps # Generate response variable
  data <- data.frame(Y = response, predictors) # Combine predictors and response into a data frame
  return(data) # Return the data frame
}

# Generate a dataset for the paper
data <- generator(n, p, rho) # Generate synthetic data
rm(generator, rho, n, p); gc() # Remove unnecessary objects and run the garbage collector

p <- 500 # Number of predictors
B <- 500 # Number of subsamples

x <- as.matrix(data[, 2:ncol(data)]) # Predictors
y <- data$Y # Response

cv_lasso <- cv.glmnet(x, y, nfolds = 10, alpha = 1) # Fit LASSO model with 10-fold CV

candidate_set <- cv_lasso$lambda # Candidate set of lambda values

S_list <- vector("list", length(candidate_set)) # Initialize a list to store selection matrix for each lambda
names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries

# Stability Selection for each lambda in candidate_set
for (lambda_idx in seq_along(candidate_set)) {
  
  lambda <- candidate_set[lambda_idx]  # Current lambda value
  S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
  colnames(S) <- colnames(x) # Set column names of S to predictor names
  
  for (i in 1:B) { 
    # Sub-sample the data (half of the original data without replacement)
    model_data <- data[sample(1:nrow(data), nrow(data) / 2, replace = FALSE), ]
    
    # Prepare the response and predictors
    x_sub <- as.matrix(model_data[, -1]) # Exclude the response variable
    y_sub <- model_data$Y # Response variable
    
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

stable_values <- which(stab_values > 0.75) # Index of stable lambda values
lambda_stable <- min(candidate_set[stable_values]) # Minimum stable lambda value
index_of_lambda_stable <- which(candidate_set == lambda_stable) # Index of lambda_stable


stability <- data.frame() # Initialize a data frame to store stability values
Stable_S <- S_list[[index_of_lambda_stable]] # Stable selection matrix for lambda_stable
for (k in 2:nrow(Stable_S)){ # loop through subsamples results
  output <- getStability(Stable_S[1:k,]) # Compute stability values
  stability <- rbind(stability, data.frame(k, output$stability, output$variance, output$lower, output$upper)) # Append stability values to the data frame
}
colnames(stability) <- c('Iteration', 'Stability', 'Variance', 'Lower', 'Upper') # Set column names of the data frame

ggplot(stability, aes(x = Iteration, y = Stability)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = 'blue', alpha = 0.7) + # Add ribbon for confidence interval
  labs(title = TeX('Stability of Stability Selection ($\\rho = 0.2$, $\\lambda = \\lambda_{stable}$)'), 
       x = 'Iteration (sub-sample)', y = TeX('Stability ($\\hat{\\Phi}$)'))+
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),       # Title text size
    axis.title.x = element_text(size = 18),     # X-axis label size
    axis.title.y = element_text(size = 18),     # Y-axis label size
    axis.text.x = element_text(size = 16),      # X-axis tick text size
    axis.text.y = element_text(size = 16)       # Y-axis tick text size
  )
