# This script generates Figure 8 in the paper.
# The script generates synthetic data with 500 predictors, 50 observations, and 
#a constant correlation coefficient of 0.2 followed by stable stability selection with Stable Lasso
# To reproduce the figure reported in the paper, you need to set the number of
# subsamples B = 500. For faster computation, you may reduce the number of subsamples.

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
n <- 50 # Number of samples
p <- 500 # Number of predictors
rho <- 0.2 # correlation coefficient

# Newton's method for optimization
newton <- function(f, fp, x, tol = 1e-3, m = 100) {
  iter <- 0
  oldx <- x
  x <- oldx + 10 * tol
  
  while (abs(x - oldx) > tol) {
    iter <- iter + 1
    if (iter > m) {
      warning("Maximum iterations (m = ", m,
              ") reached: returning last estimate", call. = FALSE)
      break
    }
    oldx <- x
    x <- x - f(x) / fp(x)
  }
  
  return(x)
}

# AirHOLP (2025) from "https://github.com/Logic314/Air-HOLP"
AirHOLP <- function(X, y, Threshold = min(dim(X)[2]-1, dim(X)[1]/log(dim(X)[1])),
                    r0 = 10, adapt = TRUE, iter = 10, Lambda, U, XU) {
  # Arguments:-
  # X: matrix of features (Matrix)
  # y: response vector (Vector)
  # Threshold: screening threshold (Integer)
  # r0: initial penalties (Vector)
  # adapt: if >= 1 adaptive penalty will be used (Binary)
  # iter: maximum number of iteration for adaptive penalty selection (Integer)
  # Lambda: eigenvalues of XXT, if missing the function will compute it (Vector)
  # U: eigenvectors of XXT, if missing the function will compute it (Matrix)
  # XU: X transpose times U, if missing the function will compute it (Matrix)
  
  # Output:-
  # index_r: ranking of features by Air-HOLP (Matrix)
  # index_r0: ranking of features by Ridge-HOLP (Matrix)
  # Beta_r: regression coefficients of Air-HOLP (Matrix)
  # Beta_r0: regression coefficients of Ridge-HOLP (Matrix)
  # r: selected penalty parameters by Air-HOLP (Vector)
  # iter_last: number of iterations used in Air-HOLP (Vector)
  
  n <- dim(X)[1] # sample size
  p <- dim(X)[2] # number of features
  q <- length(r0) # number of penalty parameters
  iter_temp2 <- 0*(1:q) # used for calculating iter_last
  iter_temp1 <- iter_temp2 - 1 # used for calculating iter_last
  
  # Standardizing X and y:
  X <- X - matrix(rep(colMeans(X),each = n),n,p)
  X <- X/matrix(rep(sqrt(colMeans(X^2)),each = n),n,p)
  y <- (y - mean(y))/sd(y)
  
  if(adapt){
    # Main computations:
    if(missing(Lambda)|missing(U)){
      XXT <- tcrossprod(X)
      eXXT <- eigen(XXT)
      Lambda <- eXXT$values
      U <- eXXT$vectors
    }
    if(missing(XU)){
      XU <- crossprod(X,U)
    }
    Dn <- diag(Lambda)
    UTy <- crossprod(U,y)
    yUD2UTy <- UTy^2*(Lambda^2)
    
    # Penalty selection:
    r_max <- 1000*sqrt(n) # maximum penalty
    max.iter <- 30 # maximum number of iterations for Newtons method
    index_r <- matrix(1:(p*q), nrow = p, ncol = q)
    index_r0 <- index_r
    Beta_r <- index_r
    Beta_r0 <- index_r
    r <- r0
    r_temp <- r0
    for (j in 1:iter) {
      for (i in 1:q) {
        # Initial screening:
        Beta_temp <- XU%*%((Lambda+r[i])^(-1)*UTy)
        index_temp <- match(1:p,rank(-abs(Beta_temp), na.last = NA,
                                     ties.method = c("random")))
        Xs <- X[,index_temp[1:Threshold]] # Screened features
        if(j<2) {
          Beta_r0[,i] <- Beta_temp
          index_r0[,i] <- rank(-abs(Beta_temp), na.last = NA,
                               ties.method = c("random"))
        }
        
        # Estimating the expected response:
        ys <- Xs%*%(solve(crossprod(Xs) +
                            diag(Threshold)*10^-12)%*%crossprod(Xs,y))
        
        # MSE functions:
        ysUDUTy <- t(crossprod(ys,U)*Lambda)*UTy
        Z <- function(lam) { # The function we minimize
          t((Lambda+lam)^-2)%*%yUD2UTy - 2*t((Lambda+lam)^-1)%*%ysUDUTy
        }
        Z1 <- function(lam) { # First derivative
          -2*t((Lambda+lam)^-3)%*%yUD2UTy + 2*t((Lambda+lam)^-2)%*%ysUDUTy
        }
        Z2 <- function(lam) { # Second derivative
          6*t((Lambda+lam)^-4)%*%yUD2UTy - 4*t((Lambda+lam)^-3)%*%ysUDUTy
        }
        
        # MSE minimization:
        sol <- newton(Z1, Z2, 0.0001, tol = 0.001, m = max.iter)
        r[i] <- sol
        if(r[i] > r_max) {r[i] <- r_max}
        if(r[i] < 0.0001) {r[i] <- 0.0001}
        if(Z(r_max) < Z(r[i])) {r[i] <- r_max} # Checking boundaries
        if(Z(0.0001) < Z(r[i])) {r[i] <- 0.0001}
        
        # Feature screening:
        Beta_r[,i] <- XU%*%((Lambda+r[i])^(-1)*UTy)
        index_r[,i] <- rank(-abs(Beta_r[,i]), na.last = NA,
                            ties.method = c("random"))
        
        # Calculations for the number of iterations:
        if(abs(r[i] - r_temp[i]) < 0.01*r[i]){ # Checking relative error
          iter_temp1[i] <- j
          iter_temp2[i] <- iter_temp2[i] + 1
        }
      }
      if(sum(abs(r - r_temp) < 0.01*r) == q){ # Checking relative error
        break
      }
      r_temp <- r
    }
    iter_last <- iter_temp1 - iter_temp2 + 1 # Number of iterations
    AirHOLP <- list(index_r = index_r, index_r0 = index_r0, Beta_r = Beta_r,
                    Beta_r0 = Beta_r0, r = r, iter_last = iter_last)
  } else{
    if(q < 2) {
      # Feature screening:
      if(missing(Lambda)|missing(U)){
        Beta_r0 <- crossprod(X, solve(tcrossprod(X)+r0*diag(n),y))
      } else{
        UTy <- crossprod(U,y)
        Beta_r0 <- XU%*%((Lambda+r0)^(-1)*UTy)
      }
      index_r0 <- rank(-abs(Beta_r0), na.last = NA, ties.method = c("random"))
      AirHOLP <- list(index_r0 = index_r0, Beta_r0 = Beta_r0)
    } else{
      # Main computations:
      if(missing(Lambda)|missing(U)){
        XXT <- tcrossprod(X)
        eXXT <- eigen(XXT)
        Lambda <- eXXT$values
        U <- eXXT$vectors
      }
      if(missing(XU)){
        XU <- crossprod(X,U)
      }
      Dn <- diag(Lambda)
      UTy <- crossprod(U,y)
      
      # Feature screening:
      index_r <- matrix(1:(p*q), nrow = p, ncol = q)
      for (i in 1:q) {
        Beta_r0[,i] <- XU%*%((Lambda+r0[i])^(-1)*UTy)
        index_r0[,i] <- rank(-abs(Beta_r0[,i]), na.last = NA,
                             ties.method = c("random"))
      }
      AirHOLP <- list(index_r0 = index_r0, Beta_r0 = Beta_r0)
    }
  }
}


# Function to generate synthetic data
generator <- function(n, p, rho){
  
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  
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


Threshold <- nrow(x)/log(nrow(x))  # Screening threshold
AHOLP <- AirHOLP(x, y, Threshold = Threshold, r0 = 10, adapt = TRUE, iter = 10)
ranked_features <- AHOLP$index_r  # Ranking of features
w <- 1 - (1/ranked_features)


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
    lasso_model <- glmnet(x_sub, y_sub, alpha = 1, lambda = lambda, penalty.factor = w)
    
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

par(mgp = c(2.2, 1, 0))  # Adjust the second value to control title spacing
plot(candidate_set, stab_values, type = "l", col = "blue", lwd = 2, 
     xlab = TeX("Regularisation Value ($\\lambda$)"), 
     ylab = TeX("Stability ($\\hat{\\Phi}$)"), 
     main = TeX("Stability vs. Regularisation Value ($\\rho = 0.2$)"),
     cex.lab = 1.5, # Increase axis title size
     cex.main = 1.8, # Increase main title size
     ylim = c(0, 1)) # Plot stability values against lambda values
abline(h = 0.75, col = "red", lty = 5)  # Add a horizontal line at stability = 0.75
abline(h = 0.4, col = "red", lty = 5)  # Add a horizontal line at stability = 0.4

index_of_min <- which(candidate_set == cv_lasso$lambda.min) # Index of lambda.min
points(candidate_set[index_of_min], stab_values[index_of_min], 
       col = "red", pch = 19, cex = 2) # Show by red dot index_of_min and stab_values[index_of_min] on the plot
text(candidate_set[index_of_min], stab_values[index_of_min], 
     "min", pos = 1, col = "red", cex = 1.5) # Add text for lambda.min

index_of_1se <- which(candidate_set == cv_lasso$lambda.1se) # Index of lambda.1se
points(candidate_set[index_of_1se], stab_values[index_of_1se], 
       col = "red", pch = 19, cex = 2) # Show by red dot index_of_1se and stab_values[index_of_1se] on the plot
text(candidate_set[index_of_1se], stab_values[index_of_1se], "1se", 
     pos = 1, col = "red", cex = 1.5) # Add text for lambda.1se

stable_values <- which(stab_values > 0.75) # Index of stable lambda values
lambda_stable <- min(candidate_set[stable_values]) # Minimum stable lambda value
index_of_lambda_stable <- which(candidate_set == lambda_stable) # Index of lambda_stable
points(candidate_set[index_of_lambda_stable], stab_values[index_of_lambda_stable], 
       col = "red", pch = 19, cex = 2) # Show by red dot index_of_lambda_stable and stab_values[index_of_lambda_stable] on the plot
text(candidate_set[index_of_lambda_stable], stab_values[index_of_lambda_stable], 
     "stable", pos = 1, col = "red", cex = 1.5) # Add text for lambda_stable



stable_values <- which(stab_values > 0.75) # Index of stable lambda values
lambda_stable <- min(candidate_set[stable_values]) # Minimum stable lambda value
index_of_lambda_stable <- which(candidate_set == lambda_stable) # Index of lambda_stable
candidate_set[index_of_lambda_stable] # Display the stable lambda value
m <- colMeans(S_list[[index_of_lambda_stable]]) # Compute selection frequencies for stable lambda
m[m > 0.6] # Display selected variables with selection probability > 0.6

m[which.max(m[m <= 0.6])] # Display the maximum selection frequency among variables with selection frequency <= 0.6
