###########################################################################
#### Simulation study ######################################################
###########################################################################

# This was repeated using random hot deck procedure, and distance hot deck 
# procedure with proxy of Cramer's V of 0.4 and 0.6

### Load libraries ###
library(tidyverse)
library(StatMatch)
library(doParallel)
library(parallel)

### Read in data ###
# data ready for sim study was saved as "sim_data.csv"
# refered to in rest of sim code as Gez_sim as was from Health data
Gez_sim <- read.csv("sim_data.csv", header = TRUE)

### Proxy variable Creation ###

# Proxy created for variable income which has 5 cateogories 1:5

# with 0.525 in diagonal, and 0.11875 elsewhere for Cramer's V = 0.397
# with 0.7 in diagonal, and 0.075 elsewhere for Cramer's V = 0.610

# Make a proxy using transition matrix

# transition matrix
trans_mat <- matrix(0.075, nrow = 5, ncol = 5)
diag(trans_mat) <- rep(0.7, times = 5)
rownames(trans_mat) <- colnames(trans_mat) <- 1:5
Gez_sim$proxy <- Gez_sim$income

# Allocated values based on transition matrix, set seed
set.seed(123)
for(i in 1:5){
  Gez_sim$proxy[Gez_sim$income == i] <- sample(1:5,
                                               size = sum(Gez_sim$income == i),
                                               prob = trans_mat[i,],
                                               replace = T)
}

# Check how related the proxy and income variables are:
# Cramer's V function
cramersV <- function(x1, x2){
  # x1 variable of interest
  # x2 variable of interest
  
  # Create the contingency table
  x1 <- as.vector(x1)
  x2 <- as.vector(x2)
  
  N <- length(x1)
  chisq <- suppressWarnings(chisq.test(x1,x2, correct = FALSE)$statistic)
  Phi <- chisq/N
  Row <- length(unique(x1))
  C <- length(unique(x2))
  
  cramersV <- sqrt(Phi/min(Row - 1, C - 1))
  
  return(cramersV)
}
# Check relation
cramersV(Gez_sim$proxy, Gez_sim$income) 


### Load in the SM method functions ###
source(paste0(scriptpath, "Code_for_the_SM_methods.R"))

### Load in quality measure ###
source(paste0(scriptpath, "Quality_measure.R"))

### Simulation Study itself ###

# Set variables names
X <- c("sex", "age", "edu", "eth", "proxy")
Y <- "income"
Z <- "health"

# Set sizes for A and B
# Calculated by the size of this sample, from Dutch pop at the time
# 421226/16979120 approx 0.0248
# 0.0248 * 421226 = 10450, so using 10500
n_A <- n_B <- 10500

# Let Y and Z be the target variables and X be matching variables

# The contingency table over the whole population:

t <- table(Gez_sim[,Y], Gez_sim[,Z])/nrow(Gez_sim)

# Make sure matching variables categorical
Gez_sim$sex <- as.character(Gez_sim$sex)
Gez_sim$age <- as.character(Gez_sim$age)
Gez_sim$edu <- as.character(Gez_sim$edu)
Gez_sim$eth <- as.character(Gez_sim$eth)
Gez_sim$proxy <- as.character(Gez_sim$proxy)

# ident is identifying variable in data

# Create lists for estimates and quality_estimates
sample_ct_A <- list()
sample_ct_B <- list()
est <- list()
quality_est_var <- list()
quality_est_rel_sd <- list()
quality_est_bias <- list()
quality_est_rel_bias <- list()
quality_est_contin <- list()
quality_est <- list()

# Code for palatalizing

# Function to put parralised data in correct format
comb <- function(x, ...){
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Set seed for sim study
set.seed(885530)

myCluster <- parallel::makeCluster(detectCores(), type = "PSOCK",
                                   output = "")
registerDoParallel(myCluster)
clusterExport(myCluster, varlist = c("n_A", "n_B", 
                                     "X", 'Y', 'Z',
                                     'dist_hotdeck', "hot_deck",
                                     'quality_SM',
                                     'Gez_sim'))

# Sim study, set number of iterations, S, with i 
simstudy <- foreach(i=1:100,
                    .packages = c("tidyverse"),
                    .combine = "comb",
                    .multicombine = TRUE,
                    .init = list(list(), list(), list(), list(),
                                 list(), list(), list(), list()))%dopar%{
                                   
                                   
                                   # Repeat iterations S times
                                   # Draw samples A_s and B_s from the complete data without replacement,
                                   # of fixed size n_A and n_B
                                   A <- Gez_sim[sample(nrow(Gez_sim), size = n_A, replace = FALSE),]
                                   B <- Gez_sim[sample(nrow(Gez_sim), size = n_B, replace = FALSE),]
                                   
                                   # Store true contingency tables of the samples, which will be matched
                                   sample_ct_A <- table(A[,Y], A[,Z])/n_A
                                   sample_ct_B <- table(B[,Y], B[,Z])/n_B
                                   
                                   # Create dataset C, from the overlap of datasets A and B, 
                                   # based on identifying variable
                                   C <- inner_join(A, B, by = c("ident", X, Y, Z))
                                   
                                   # Remove observations for Z in A
                                   A <- A[,c(Y, X, "ident")]
                                   # Remove observations for Y in B
                                   B <- B[,c(Z, X, "ident")]
                                   
                                   # Perform Statistical matching on A and B
                                   SM <- SM_method(A = A, B = B, X = X, Y = Y, Z = Z)
                                   # where SM_method is the desired SM method function used in this sim study
                                   # In my project with hot_deck for random hot deck procedure
                                   # or dist_hotdeck for distance hot deck procedure
                                   
                                   # Estimate contingency table from the matched data,
                                   # which will be used to estimate bias and variance
                                   est <- table(SM[,Y], SM[,Z])/nrow(SM)
                                   
                                   # Apply proposed bootstrap procedure on samples A and B
                                   
                                   # Make variables for A and B without values in C
                                   # Find values in A and B which are in C
                                   CinA <- which(A$ident %in% C$ident) 
                                   CinB <- which(B$ident %in% C$ident) 
                                   
                                   # Perform proposed quality procedure, removing values in C from A and B
                                   # Perform with just matching and target variables
                                   # R = 200
                                   quality_est <- quality_SM(A = A[-CinA, c(Y,X)], 
                                                             B = B[-CinB, c(Z,X)], 
                                                             C = C[,c(Z,Y,X)],
                                                             SM = SM_method, Y = Y, Z = Z, X = X,
                                                             R = 200, aux = 1)
                                   # where SM_method is the desired SM method function used in this sim study
                                   # In my project with hot_deck for random hot deck procedure
                                   # or dist_hotdeck for distance hot deck procedure
                                   # Make sure they are the same as above
                                   
                                   # Store results of the quality procedure so they can be return separately
                                   quality_est_bias <- quality_est[[1]]
                                   quality_est_rel_bias <- quality_est[[2]]
                                   quality_est_var <- quality_est[[3]]
                                   quality_est_rel_sd <- quality_est[[4]]
                                   quality_est_contin <- quality_est[[5]]
                                   
                                   # return all relevant values
                                   return(list(sample_ct_A,
                                               sample_ct_B,
                                               est,
                                               quality_est_bias,
                                               quality_est_rel_bias,
                                               quality_est_var,
                                               quality_est_rel_sd,
                                               quality_est_contin))
                                   
                                 }
stopCluster(myCluster)

# Set S, the number of sim iterations
S <- 100

# Mean of simulation estimates
mean_est <- (1/S)*Reduce('+', est)

# Estimate for true bias 
true_bias_list <- lapply(est, FUN = function(x){x - t})
true_bias_red <- Reduce('+', true_bias_list)
true_bias <- (1/S)*true_bias_red

# Estimate for true relative bias
true_rel_bias <- (mean_est - t)/t

# Estimate for true variance 
true_var_list <- lapply(est, FUN = function(x){(x - mean_est)^2})
true_var_red <- Reduce('+', true_var_list)
true_variance <- (1/(S-1))*true_var_red

# True relative standard deviation estimate
true_rel_sd <- sqrt(true_variance)/mean_est

#### Compare bootrap and true estimates ### 

### Averages of the bootstrap estimates ###
# Bootstrap bias estimate average
bias_avg <- (1/S)*Reduce('+', quality_est_bias)

# Bootstrap relative bias estimate average 
rel_bias_avg <- (1/S)*Reduce('+', quality_est_rel_bias)

# Bootstrap variance estimate average
var_avg <- (1/S)*Reduce('+', quality_est_var)

# Bootstrap relative sd estimate average
rel_sd_avg <- (1/S)*Reduce('+', quality_est_rel_sd)


### Absolute differences between bootstrap and true estimates ###
# Absolute differences in bootstrapped estimates and the truth for bias
abs_bias_list <- lapply(quality_est_bias, 
                        FUN = function(x){abs(x - true_bias)})
abs_bias <- (1/S)*Reduce('+', abs_bias_list)

# Absolute differences in bootstrapped estimates and the truth for relative
# bias
abs_bias_rel_list <- lapply(quality_est_rel_bias, 
                            FUN = function(x){abs(x - true_rel_bias)})
abs_bias_rel <- (1/S)*Reduce('+', abs_bias_rel_list)

# Absolute differences in bootstrapped estimates and the truth for relative
# sd
abs_sd_rel_list <- lapply(quality_est_rel_sd, 
                          FUN = function(x){abs(x - true_rel_sd)})
abs_sd_rel <- (1/S)*Reduce('+', abs_sd_rel_list)

# Absolute differences in bootstrapped estimates and the truth for variance
abs_var_list <- lapply(quality_est_var, 
                       FUN = function(x){abs(x - true_variance)})
abs_var <- (1/S)*Reduce('+', abs_var_list)

### The bootstrap estimates divided by the true estiamtes ###
# Divided by truth for bias
bias_div_list <- lapply(quality_est_bias, 
                        FUN = function(x){x/true_bias})
bias_div <- (1/S)*Reduce('+', bias_div_list)

# Divided by truth for variance
var_div_list <- lapply(quality_est_var, 
                       FUN = function(x){x/true_variance})
var_div <- (1/S)*Reduce('+', var_div_list)

### Estimation error calculation
# True estimation error

# Calculate the total contingency tables for samples A and B
sample_ct <- list()
for(i in 1:100){
  sample_ct[[i]] <- ((sample_ct_B[[i]]*n_B)+(sample_ct_A[[i]]*n_A))/
    (n_A + n_B)
}
true_est_err_list <- Map("-", est, sample_ct)
true_est_err <- (1/S)*Reduce('+', true_est_err_list)

# difference between true estimation error and true bias
true_est_err - true_bias

# True relative estimation error
true_rel_err_list <- Map("/", true_est_err_list, sample_ct)
true_rel_err <- (1/S)*Reduce('+', true_rel_err_list)

# difference between relative estimation error and 
# true relative bias
true_rel_err - true_rel_bias

