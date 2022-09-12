####################################################
#### Codes for the Statistical Matching Methods ####
####################################################


#################################
### Random Hot deck Procedure ###
#################################

# Load Library
# library(tidyverse)

### Random hot deck procedure function ###

hot_deck <- function(A, B, X, Y, Z){
  # A - dataset for sample A
  # B - dataset for sample B
  # Where dataset B will be matched with values from dataset A
  # (A is donor and B is receiving dataset)
  # X - the names in a categorical list for matching variables 
  # Matching variables assumed categorical
  # these must be labeled the same for both dataset_A and dataset_B
  # Y - the variable name for one target variable in dataset A
  # Z - the variable name for the other target variable in dataset B
  
  # require
  require(dplyr)
  
  # Check that the number of observations in is not A is less than B
  if(nrow(A) < nrow(B)){
    cat("Warning: the number of observations in A is less than the number
        of oberservations in B", fill = TRUE)
  }
  
  # Make sure all matching variables are categorical 
  A[,X] <- sapply(A[,X], as.character)
  B[,X] <- sapply(B[,X], as.character)
  # Make sure it is a dataframe
  A <- as.data.frame(A)
  B <- as.data.frame(B)
  
  # Create dataframe for the output
  B_matched <- B
  B_matched$Y <- rep(NA, times = nrow(B_matched))
  colnames(B_matched)[colnames(B_matched) == 'Y'] <- Y
  
  # Variable for 'donor' unit location in A
  B_matched$matched_donor <- rep(NA, times = nrow(B_matched))
  
  # Groups within B dataframe
  # Add columns to dataframe of B, indicating a group number based on 
  # the matching variables
  B <- B %>% 
    group_by(across(all_of(X))) %>% 
    mutate(group = cur_group_id())
  options(dplyr.summarise.inform = FALSE)
  # Create a dataframe of the groups of X variables that are indicative 
  # of the group numbers in previous column creation
  groupsB <- B %>% group_by(across(all_of(X))) %>%
    summarise() %>%
    rowwise() 
  groupsB <- as.data.frame(groupsB)
  
  # A table of all combinations of variables, that can be removed
  comb <- plyr::ldply(1:length(X), function(x)t(combn(1:length(X), x)))
  comb[is.na(comb)] <- 0
  # value which denotes which combination of variables will be removed
  remove <- numeric(1)
  
  
  # Assign values of A to group B, based on groups.
  # Start by see if each group in B has at least one value
  # in A from same group - if not remove variables until group match
  for(i in 1:nrow(groupsB)){
    # Starting value so no variables removed initially, 
    # and records which groups will be
    var <- 0
    # Variable to log if a group match has been found in A
    pairing <- numeric(1)
    # Set grouping, if becomes 2, then no groups used for sampling
    grouping <- 1
    # A value of recording the names of X in useage
    using_var <- X
    # A version of A dataset with only the X variables in use
    A_use <- A
    
    
    while(pairing < 1){
      # If combination wasn't found initially var will be set above 0
      if(var != 0){ 
        # Values to be removed in this round
        remove <- unlist(comb[var,])
        # The variables to be used in this round
        using_var <- X[-remove]
        # Create new version of A_use
        A_use <- A[,c(Y, using_var)]
      }
      
      # The group category of B trying to be matched in this for loop iteration
      # Is restricted when using less variables
      group <- groupsB[i,using_var]
      
      # Make variable in A_use of the groups in A to cross-reference with B
      A_use <- A_use %>% 
        group_by(across(all_of(using_var))) %>% 
        mutate(group = cur_group_id())
      
      # Data frame which is reference to te grouping in A
      groupsA <- A_use %>% group_by(across(all_of(using_var))) %>%
        summarise() %>%
        rowwise()
      groupsA <- as.data.frame(groupsA)
      
      # See if B group used in this round is in the A groups
      if(sum(duplicated(rbind(group, groupsA))) < 1){
        # Here B group does not have match in A, so increase var
        # to try and find match
        var <- var + 1
        # Check if all combination of matching variables have been tried
        # If so set to use random matching across dataset
        if(var > sum(sapply(1:(length(X)-1),
                            function(x)choose(length(X),x)))){
          # No combination possible with matching variables so
          # will use random hot deck
          pairing <- 1
          grouping <- 2
        }
      }else{# Here a match is found
        pairing <- 1}
      
    }
    
    # observations of group i in B
    Bi <- which(B$group == i)
    
    if(grouping == 1){
      # Which group in A_use is matching group i in B
      Agroup <- which(duplicated(rbind(group, groupsA))) - 1
      
      # Find observations in A that can be matched to this group in B
      Amatch <- A_use
      Amatch$row_num <- 1:nrow(Amatch)
      Amatch <- Amatch[which(Amatch$group == Agroup),]
      
      # Match these values to group i values in B
      for(j in Bi){
        match <- sample(1:nrow(Amatch), size = 1)
        B_matched[j,'matched_donor'] <- Amatch[match, "row_num"]
        B_matched[j, Y] <- Amatch[match, Y]
      }
    }
    
    if(grouping == 2){ # where no combination was found
      # so random matching is used 
      
      # Possible matching dataset
      Amatch <- A_use
      Amatch$row_num <- 1:nrow(Amatch)
      
      for(j in Bi){
        # Random hotdeck from rest of datasets for group i values in B
        match <- sample(1:nrow(A), size = 1)
        B_matched[j,'matched_donor'] <- Amatch[match, "row_num"]
        B_matched[j, Y] <- Amatch[match, Y]
      }
    }
  }
  
  
  # Return the matched dataset with matching and target variables
  return(B_matched)
  
}

### Distance hot-deck procedure using StatMatch ###

dist_hotdeck <- function(A, B, X, Y, Z){
  # A - dataset for sample A
  # B - dataset for sample B
  # Where dataset B will be matched with values from dataset A
  # (A is donor and B is receiving dataset)
  # X - the names in a categorical list for matching variables 
  # Matching variables assumed categorical
  # these must be labeled the same for both dataset_A and dataset_B
  # Y - the variable name for one target variable in dataset A
  # Z - the variable name for the other target variable in dataset B
  
  require(StatMatch)
  
  # Datasets for use in the matching 
  A_use <- A[,c(X, Y)]
  B_use <- B[,c(X, Z)]
  
  # The function used in statmatch to match datasets, using Gower distance
  # dataset B in receiveng and A is the donor, X are matching variables
  stat_match <- NND.hotdeck(data.rec = B_use, data.don = A_use, 
                            match.vars = X,
                            dist.fun = "Gower")
  # Creates it into usable dataframe
  data <- create.fused(data.rec=B_use, data.don=A_use,
                       mtc.ids=stat_match$mtc.ids, z.vars=Y)
  
  # return dataframe
  return(data)
}



