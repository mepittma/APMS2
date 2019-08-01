####FUNCTIONS####
remove_low_replicate <- function(data, minReplicates=2, missing=NA, replaceWith=0){
  #replaces with 0 every value that is nonzero in less than the required number of replicates
  
  samples <- colnames(data)
  variant <- gsub("_.$", "", samples)
  replicate = c()
  idx = 1
  for(i in range(1,length(unique(variant)))){
    n = sum(variant==variant[idx])
    replicate = c(replicate, rep(i,n))
    idx = idx + n
  }
  
  groups <- split(1:ncol(data), f = as.vector(replicate))
  for(x in groups){
    #print(head(data[x]))
    if(is.na(missing)){
      nonZeros <- apply(data[x], MARGIN = 1, FUN = function(x){sum(!is.na(x))})}
    else if(missing==0){nonZeros <- apply(data[x], MARGIN = 1, FUN = function(x){sum(x!=0)})}
    #print(head(nonZeros))
    data[x][nonZeros<minReplicates,] <- replaceWith
    #print(head(data[x]))
  }
  #print(head(data))
  return(data)
}

impute_data = function(df, width = 0.5, downshift = 3, missingVal=NA, byColumn=FALSE) {
  # df = data frame containing filtered 
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  dfOriginal <- as.matrix(df)
  dfNames = names(df)
  imputeNames = paste0(names(df), "_impute")
  
  # Create new columns indicating whether the values are imputed
  df[imputeNames] = lapply(names(df), function(x) df[, x] %in% c(missingVal))
  
  # Imputation
  set.seed(100)
  # print(head(df))
  # print(table(df==0))
  if(byColumn==FALSE){
    print("IMPUTING values by total matrix")
    tempSD <- width * sd(as.matrix(dfOriginal[!dfOriginal %in% c(NA, 0)], na.rm = TRUE))
    tempMean = mean(as.matrix(dfOriginal[!dfOriginal %in% c(NA, 0)], na.rm = TRUE)) - 
      downshift * sd(as.matrix(dfOriginal[!dfOriginal %in% c(NA, 0)], na.rm = TRUE))   
    # shift mean of imputed values
    nMissing = sum(dfOriginal %in% c(missingVal))
    print(tempSD)
    print(tempMean)
    print(rnorm(nMissing, mean = tempMean, sd = tempSD))
    if(nMissing!=0){dfOriginal[dfOriginal %in% c(missingVal)] <-  rnorm(nMissing, mean = tempMean, sd = tempSD)}
    return(dfOriginal)
    stop()
  }
  if(byColumn==TRUE){
    print("IMPUTING values by column")
    df[dfNames] = lapply(dfNames,
                         function(x) {
                           #print(x)
                           temp <- df[[x]]
                           #print(temp)
                           
                           tempSD <- width * sd(temp[!temp %in% c(NA, 0)], na.rm = TRUE)   # shrink sd width
                           tempMean = mean(temp[!temp %in% c(NA, 0)], na.rm = TRUE) - 
                             downshift * sd(temp[!temp %in% c(NA, 0)], na.rm = TRUE)   
                           # shift mean of imputed values
                           nMissing = sum(temp %in% c(missingVal))
                           
                           #print(tempSD)
                           #print(tempMean)
                           # print(nMissing)
                           # print(rnorm(nMissing, mean = tempMean, sd = tempSD))
                           # print(temp[temp %in% c(missingVal)])
                           if(nMissing != 0){temp[temp %in% c(missingVal)] <-  rnorm(nMissing, mean = tempMean, sd = tempSD)}
                           # print(temp)
                           return(temp)
                         })
    print("Printing Plot")
    #Plots distribution of log values for each sample (column)
    #panes <- c(round(ncol(df)/2/4+0.1), 4)
    #dfLong <- cbind(melt(as.matrix(df[dfNames])), melt(as.matrix(df[imputeNames]))['value'])
    #names(dfLong) <- c("Protein", "Sample", "Intensity", "Imputed")
    #dfLong$rep <- gsub(".*_", "", dfLong$Sample)
    #dfLong$bait <- gsub("_.$", "", dfLong$Sample)
    #ggplot(data = dfLong, aes(x=Intensity, fill=Imputed))
    #+geom_density(alpha=0.3)+facet_grid(bait ~ rep)
    #+ggtitle(label="Distribution of Imputed vs Empirical values")
    #print(head(df[!grepl("impute", names(df))]))
    return(df[!grepl("impute", names(df))])
  }
}
