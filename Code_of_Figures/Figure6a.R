# This script generates Figure 6a in the paper.
# The script generates synthetic data with 500 predictors, 
# 50 observations for training, and 25 observations for testing. 
# correlation coefficient is set to 0.2.
# To reproduce the figure reported in the paper, 
#you need to set the number of
# subsamples B = 500. For faster computation, 
#you may reduce the number of subsamples.

options(warn=-1) # Turn off warnings
if (!requireNamespace("MASS")) {install.packages("MASS")} # install package if not already installed
if (!requireNamespace("glmnet")) {install.packages("glmnet")} # install package if not already installed
if (!requireNamespace("latex2exp")) {install.packages("latex2exp")} # install package if not already installed

# Load necessary library
library(MASS) # load the MASS package for mvrnorm
library(glmnet) # load the glmnet package for LASSO
library(latex2exp) # load the latex2exp package for TeX

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
#  [1] latex2exp_0.9.6 glmnet_4.1-8    Matrix_1.6-1.1  MASS_7.3-60    
#loaded via a namespace (and not attached):
#  [1] codetools_0.2-20  shape_1.4.6.1     lattice_0.22-6    magrittr_2.0.3   
#[5] splines_4.3.2     glue_1.8.0        iterators_1.0.14  stringr_1.5.1    
#[9] lifecycle_1.0.4   cli_3.6.3         foreach_1.5.2     grid_4.3.2       
#[13] compiler_4.3.2    rstudioapi_0.16.0 tools_4.3.2       Rcpp_1.0.13      
#[17] survival_3.7-0    rlang_1.1.4       stringi_1.8.4 

set.seed(26) # Set seed for reproducibility
n <- 50 # Number of training samples
p <- 500 # Number of predictors
rho <- 0.2 # correlation coefficient
B <- 500 # Number of subsamples

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

train <- generator(n, p, rho) # Generate synthetic training data
x <- as.matrix(train[, 2:ncol(train)]) # Predictors
y <- train$Y # Response

cv_lasso <- cv.glmnet(x, y, nfolds = 10, alpha = 1) # Fit LASSO model with 10-fold CV

candidate_set <- cv_lasso$lambda # Candidate set of lambda values

S_list <- vector("list", length(candidate_set)) # Initialize a list to store selection matrix for each lambda
names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries

set.seed(27) # Set seed for reproducibility
n <- 25 # Number of testing samples
test_data <- generator(n, p, rho) # Generate test data

x_test <- as.matrix(test_data[, -1]) # Test predictors
y_test <- test_data$Y # Test response

MSE <- c() # Initialize a vector to store MSE values
set.seed(26) # Set seed for reproducibility
# Stability Selection for each lambda in candidate_set
for (lambda_idx in seq_along(candidate_set)) {
  lambda <- candidate_set[lambda_idx]  # Current lambda value
  S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
  colnames(S) <- colnames(x) # Set column names of S to predictor names
  mse_value <- c() # Initialize a vector to store MSE values
  for (i in 1:B) {
    # Sub-sample the training data (half of the original data without replacement)
    model_data <- train[sample(1:nrow(train), nrow(train) / 2, replace = FALSE), ]
    
    x_sub <- as.matrix(model_data[, -1]) # Exclude the response variable
    y_sub <- model_data$Y # Response variable
    
    # Fit the LASSO model with the current lambda
    lasso_model <- glmnet(x_sub, y_sub, alpha = 1, lambda = lambda)
    
    y_pred <- predict(lasso_model, s = lambda, newx = x_test) # Predict the response for the test data
    mse <- mean((y_test - y_pred)^2) # Calculate the mean squared error
    mse_value <- c(mse_value, mse) # Store the MSE value
    # Extract significant predictors (ignoring the intercept, hence [-1])
    significant_predictors <- ifelse(coef(lasso_model) != 0, 1, 0)[-1]
    
    S[i, ] <- significant_predictors # Store the significant predictors in matrix S
  }
  MSE[lambda_idx] <- mean(mse_value) # Store the average MSE for the current lambda
  S_list[[lambda_idx]] <- S # Store the selection matrix for the current lambda
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


stability_results <- lapply(S_list, function(S) {
  getStability(S) # Compute stability values for each lambda
})

stab_values <- c() # Initialize a vector to store stability values
for (i in 1:length(candidate_set)) {
  temp <- stability_results[[paste0("lambda_", i)]]$stability # Extract stability values
  stab_values <- c(stab_values, temp) # Store stability values in the vector
}


# Pareto 
# Calculate accuracy as the negative of MSE
accuracy <- -MSE

pareto_criteria <- function(stab_values, accuracy) { # Function to find Pareto optimal candidates
  is_pareto <- rep(TRUE, length(stab_values))
  for (i in 1:length(stab_values)) {
    for (j in 1:length(stab_values)) {
      if (stab_values[j] > stab_values[i] & accuracy[j] > accuracy[i]) {
        is_pareto[i] <- FALSE
        break
      }
    }
  }
  return(which(is_pareto))
}

par(mgp = c(2.2, 1, 0))  # Adjust the second value to control title spacing
plot(candidate_set, stab_values, type = "l", col = "blue", lwd = 2, 
     xlab = TeX("Regularisation Value ($\\lambda$)"), 
     ylab = TeX("Stability ($\\hat{\\Phi}$)"), 
     main = TeX("Stability vs. Regularisation Value vs. MSE ($\\rho = 0.2$)"),
     cex.lab = 1.5, # Increase axis title size
     cex.main = 1.8, # Increase main title size
     ylim = c(0, 1))
abline(h = 0.75, col = "red", lty = 5)  # Add a horizontal line at stability = 0.75
abline(h = 0.4, col = "red", lty = 5)  # Add a horizontal line at stability = 0.4
index_of_min <- which(candidate_set == cv_lasso$lambda.min) # Index of lambda.min
points(candidate_set[index_of_min], stab_values[index_of_min], col = "red", pch = 19, cex = 2) # Show by red dot index_of_min and stab_values[index_of_min] on the plot
text(candidate_set[index_of_min], stab_values[index_of_min], "min", pos = 1, col = "red", cex = 1.5) # Add text lambda.min

index_of_1se <- which(candidate_set == cv_lasso$lambda.1se) # Index of lambda.1se
points(candidate_set[index_of_1se], stab_values[index_of_1se], col = "red", pch = 19, cex = 2) # Show by red dot index_of_1se and stab_values[index_of_1se] on the plot
text(candidate_set[index_of_1se], stab_values[index_of_1se], "1se", pos = 1, col = "red", cex = 1.5) # Add text 1se

stable_values <- which(stab_values > 0.75) # Index of stability > 0.75
lambda_stable <- min(candidate_set[stable_values]) # min lambda for stability > 0.75
index_of_lambda_stable <- which(candidate_set == lambda_stable) # Index of lambda_stable
points(candidate_set[index_of_lambda_stable], stab_values[index_of_lambda_stable], col = "red", pch = 19, cex = 2) # Show by red dot index_of_lambda_stable and stab_values[index_of_lambda_stable] on the plot
text(candidate_set[index_of_lambda_stable], stab_values[index_of_lambda_stable], "stable", pos = 1, col = "red", cex = 1.5) # Add text stable

pareto_indices <- pareto_criteria(stab_values, accuracy) # Find Pareto optimal candidates
#index_of_lambda_stable %in% pareto_indices
#[1] TRUE
# which means that lambda_stable is a Pareto optimal candidate
pareto_sums <- stab_values[pareto_indices] + accuracy[pareto_indices] # Calculate the sum of stability and accuracy for Pareto optimal candidates

index_of_lambda_pareto <- pareto_indices[which.max(pareto_sums)] # Index of lambda_pareto

points(candidate_set[index_of_lambda_pareto], stab_values[index_of_lambda_pareto], col = "red", pch = 19, cex = 2) # Show by red dot index_of_lambda_pareto and stab_values[index_of_lambda_pareto] on the plot
text(candidate_set[index_of_lambda_pareto], stab_values[index_of_lambda_pareto], "Pareto", pos = 1, col = "red", cex = 1.5) # Add text Pareto

par(new = TRUE)  # Add a new plot on top of the existing plot
plot(candidate_set, MSE, type = "l", col = "black", lwd = 2, axes = FALSE,
     xlab = "", ylab = "", ylim = range(MSE)) # Plot MSE values
axis(side = 4, col = "black", col.axis = "black")  # Add secondary y-axis on the right
mtext("MSE", side = 4, line = 0.5, col = "black", cex = 1.5)  # Label the secondary y-axis

points(candidate_set[index_of_min], MSE[index_of_min], col = "black", pch = 19, cex = 2) # Show by black dot index_of_min and MSE[index_of_min] on the plot
points(candidate_set[index_of_1se], MSE[index_of_1se], col = "black", pch = 19, cex = 2) # Show by black dot index_of_1se and MSE[index_of_1se] on the plot
points(candidate_set[index_of_lambda_stable], MSE[index_of_lambda_stable], col = "black", pch = 19, cex = 2) # Show by black dot index_of_lambda_stable and MSE[index_of_lambda_stable] on the plot
points(candidate_set[index_of_lambda_pareto], MSE[index_of_lambda_pareto], col = "black", pch = 19, cex = 2) # Show by black dot index_of_lambda_pareto and MSE[index_of_lambda_pareto] on the plot

abline(v = candidate_set[index_of_min], col = "red", lty = 5)  # Add a vertical line at lambda.min
abline(v = candidate_set[index_of_1se], col = "red", lty = 5)  # Add a vertical line at lambda.1se
abline(v = candidate_set[index_of_lambda_stable], col = "red", lty = 5)  # Add a vertical line at lambda_stable
abline(v = candidate_set[index_of_lambda_pareto], col = "red", lty = 5) # Add a vertical line at lambda_pareto

# Adding a legend for the two lines
legend("bottomright", legend = c("Stability", "MSE"),
       col = c("blue", "black"), 
       lty = 1, 
       lwd = 2, 
       bty = "n")  # No box around the legend