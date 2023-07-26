
library(huge)
library(tidyverse)

meanImpute = function(x)
{
  x[is.na(x)] <- mean(x,na.rm=T)
  return(x)
}

# this function removes genes that weren't measured anywhere and should be done agnostic to subtype
removeMissing = function(meth_df, thres = 0.2)
{
  propMiss = apply(meth_df,2,function(x){sum(is.na(x))/length(x)})
  skipIndices = which(propMiss > thres)
  print(paste("Number of genes with > 0.2 missing:",length(skipIndices)))
  print("Genes omitted:")
  print(names(meth_df)[skipIndices])
  meth_dfClean = meth_df[,-skipIndices]
  return(meth_dfClean)
}

betaToM = function(beta)
{
  return(log2(beta/(1-beta))) # reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
}

# this function does imputation and transformation and should be applied within subtype
cleanMethylationData = function(meth_df, npn=T, mval=F) # meth_df is a data frame of beta means, 
# rows=samples, first column=sample ids,  cols=genes
{
  meth_dfWhole = meth_df#removeMissing(meth_df,thres = 0.000000000000000000001)
  print(head(meth_df))
  print(head(meth_dfWhole))
  # skip anything w/more than 0.2 missing
  
  meth_dfComplete = meth_dfWhole
  for(i in 2:ncol(meth_dfComplete))
  {
    meth_dfComplete[,i] = meanImpute(meth_dfWhole[,i])
  }
  
  # sanity check the imputation
  summary(apply(meth_dfComplete,2,function(x){sum(is.na(x))}))
  summary(apply(meth_dfComplete[,-1],2,sd))
  
  if(mval & !npn)
  {
    transf_data = data.frame(apply(meth_dfComplete[,-1],2,betaToM))
    row.names(transf_data)= meth_dfComplete[,1]
    return(transf_data)
  }
  
  if(npn)
  {
    # do the nonparanormal transformation
    transf_data = data.frame(huge.npn(meth_dfComplete[,-1]))
    row.names(transf_data)= meth_dfComplete[,1]
    return(transf_data)
  }
  
  if(mval & npn)
  {
    # first apply M value transformation
    transf_data_m = apply(meth_dfComplete[,-1],2,betaToM)
    # next apply npn
    transf_data_m_npn = data.frame(huge.npn(transf_data_m))
    row.names(transf_data_m_npn)= meth_dfComplete[,1]
    return(transf_data_m_npn)
  }
  
  row.names(meth_dfComplete) = meth_dfComplete[,1]
  
  return(meth_dfComplete[,-1])
  
}
