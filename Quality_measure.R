######################################
########## Quality measure ###########
######################################

# Quality measure function
quality_SM <- function(A, B, C, SM, Y, Z, X, R, aux = 1){
  # A - dataset A that has obvs for target variable Y, and matching variables X
  # B - dataset B that has obvs for target variable Z, and matching variables X
  # C - dataset C that has obvs for Y, Z and X
  # C is overlap of A and B, but A and B should not have values from C in 
  # them when entered into function
  # It is assumed that target variables Y and Z are categorical
  # SM - the statistical matching method
  # Y - variable name for target variable Y
  # Z - variable name for target variable Z
  # X - variables names for matching variables
  # R - number of bootstrap samples wanted
  # This code assumes the SM method does not use dataset C explicitly,
  # and that has a function
  # in the formation of SM(A, B, X, Y, Z)
  
  # Make sure that target variables are coded as categorical
  A[,Y] <- sapply(A[,Y], as.character)
  B[,Z] <- sapply(B[,Z], as.character)
  C[,c(Y, Z)] <- sapply(C[,c(Y,Z)], as.character)
  
  # Make sure that A, B, C are dataframes
  A <- as.data.frame(A)
  B <- as.data.frame(B)
  C <- as.data.frame(C)
  
  # record the number of observations in the datasets
  na <- nrow(A)
  nb <- nrow(B)
  nc <- nrow(C)
  
  # Create contingency table from dataset C
  contin_table <- table(C[,Y], C[,Z])/nc
  
  # Check that the contingency table has all the categorical values of Y and Z
  # found in A and B, as they may not be in C
  # If they are not found, add respective rows or columns to contingency table
  # Values of Y
  if(nrow(contin_table) != length(unique(A[,Y]))){
    # Add needed rows
    rowsneed <- which(!(unique(A[,Y]) %in% rownames(contin_table)))
    if(length(rowsneed == 1)){
      rows <- t(as.matrix(replicate(ncol(contin_table),
                                    rep(0, times = length(rowsneed)))))
    }else{
      rows <- as.matrix(replicate(ncol(contin_table),
                                  rep(0, times = length(rowsneed))))
    }
    rownames(rows) <- unique(A[,Y])[rowsneed]
    contin_table <- rbind(contin_table, rows)
    # Sort in order, like table function does
    contin_table <- contin_table[sort(as.character(unique(A[,Y]))),]
  }
  #Values of Z
  if(ncol(contin_table) != length(unique(B[,Z]))){
    # Add needed rows
    colsneed <- which(!(unique(B[,Z]) %in% rownames(contin_table)))
    cols <- as.matrix(replicate(length(colsneed),
                                rep(0, times = nrow(contin_table))))
    colnames(cols) <- unique(B[,Z])[colsneed]
    contin_table <- cbind(contin_table, cols)
    # Sort in order, like table function does
    contin_table <- contin_table[,sort(as.character(unique(B[,Z])))]
  }
  
  # Create list to store the contingency tables of the bootstrap samples
  boot_contin <- list()
  
  # The bootstrapping iterations, done R times
  for(i in 1:R){
    # Draw na + nc values of cateogries of Y based on contingency table of C
    Yvals <- sample(rownames(contin_table), size = (na + nc),
                    replace = TRUE, prob = rowSums(contin_table))
    # Allocate whether they are from dataset A or C based on probability of 
    # A of (na)/(na+nc) and probability of C of (nc)/(na+nc)
    AC <- sample(c("A","C"), size = (na+nc), replace = T,
                 prob = c((na)/(na+nc), (nc)/(na+nc)))
    # Create table to show how many values of each category will come from A
    # and how many from C
    YvalsAC <- table(Yvals, AC)
    
    # Create data sets for bootstrap samples A and C
    Asamp <- NULL
    Csamp <- NULL
    
    # Assign the values from A and C to bootstrap sample
    # with each category of Y based on the table above, with replacement
    for(j in 1:nrow(YvalsAC)){
      Asamp <- rbind(A[sample(which(A[,Y] == rownames(YvalsAC)[j]), 
                              size = YvalsAC[j,1], replace = TRUE),], Asamp)
      Csamp <- rbind(C[sample(which(C[,Y] == rownames(YvalsAC)[j]), 
                              size = YvalsAC[j,2], replace = TRUE),], Csamp)
    }
    
    # Now to get bootstrap sample B
    # B will be of size db where db = nb + nc - dc
    # with dc the number of obvervations in Csamp
    dc <- nrow(Csamp)
    
    # Draw nb + nc - dc times the categories of Z based on the column sums
    # of the contingency table of C
    Zvals <- sample(colnames(contin_table), size = (nb + nc - dc),
                    replace = TRUE, prob = colSums(contin_table))
    # Create table of how many of each category will be drawn from B
    Ztots <- table(Zvals)
    
    # Create dataset for bootstrap sample B
    Bsamp <- NULL
    
    # Assign values from B to bootstrap sample based on the number of allocated 
    # in table above and the respective Z values in B
    for(j in 1:nrow(Ztots)){
      Bsamp <- rbind(B[sample(which(B[,Z] == rownames(Ztots)[j]), 
                              size = Ztots[j], replace = TRUE),], Bsamp)
    }
    
    ## Apply the SM method to the bootstrapped samples
    # including the C sample in both the A and B samples
    SM_data <- SM(A = rbind(Asamp, Csamp[,c(X, Y)]),
                  B = rbind(Bsamp, Csamp[,c(X, Z)]),
                  Y = Y, Z = Z, X = X)
    
    ## Create contingency table for matched data of bootstrap sample
    table_boot <- table(SM_data[,Y], SM_data[, Z])/nrow(SM_data)
    
    # Check it has all rows and columns as the contingency table
    # therefore having all the categories of Y and Z
    # order the contingecy table, in same order as C contingency table
    if(!all(dim(table_boot) == dim(contin_table))){
      if(nrow(table_boot) != nrow(contin_table)){
        rowsneed <- which(!(rownames(contin_table) %in% rownames(table_boot)))
        rows <- as.matrix(replicate(ncol(table_boot),
                                    rep(0, times = length(rowsneed))))
        if(length(rowsneed == 1)){
          rows <- t(as.matrix(replicate(ncol(table_boot),
                                        rep(0, times = length(rowsneed)))))
        }else{
          rows <- as.matrix(replicate(ncol(table_boot),
                                      rep(0, times = length(rowsneed))))
        }
        rownames(rows) <- rownames(contin_table)[rowsneed]
        table_boot <- rbind(table_boot, rows)
        table_boot <- table_boot[rownames(contin_table),]
      }
      if(ncol(table_boot) != ncol(contin_table)){
        colsneed <- which(!(colnames(contin_table) %in% colnames(table_boot)))
        cols <- as.matrix(replicate(length(colsneed),
                                    rep(0, times = ncol(table_boot))))
        colnames(cols) <- colnames(table_boot)[colsneed]
        table_boot <- cbind(table_boot, cols)
        table_boot <- table_boot[,colnames(contin_table)]
      }
    }
    
    # Add this iteration contingency table to boot_contin
    # boot_contin are the contingency tables from the bootstraps
    boot_contin[[i]] <- table_boot
  }
  
  # Calculate mean value of contingency tables of bootstrap samples
  means <- (1/R)*Reduce('+', boot_contin)
  
  # Bias estimation
  bias <- means - contin_table
  # relative bias estimation
  rel_bias <- (means - contin_table)/contin_table
  
  # Variance estimation
  var_list <- lapply(boot_contin, FUN = function(x){(x - means)^2})
  variance <- (1/(R-1))*Reduce('+', var_list)
  
  # Relative standard deviation estiation
  rel_sd <- sqrt(variance) / means
  
  # return all these values in a list
  boot <- list(bias, rel_bias, variance, rel_sd, boot_contin)
  
  return(boot)
}