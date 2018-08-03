sequence <- function(coverage, SNPindex, dfstrains, strainweight) {
  sequenced <- data.frame(CHR = SNPindex$CHR , POS = SNPindex$POS, depth1 = NA, ref1 = NA, alt1 = NA,depth2 = NA, ref2 = NA, alt2 = NA )
  coverage1 = rpois(n = 1, lambda = coverage)
  coverage2 = rpois(n = 1, lambda = coverage)
  for (i in 1:nrow(dfstrains)) {
    
    
    depth1 = sample(rpois(n = 1, lambda = coverage1), size =1)
    depth2 = sample(rpois(n = 1, lambda = coverage2), size =1)
    
    sequenced$depth1[i]= depth1
    sequenced$depth2[i]= depth2
    
    a<- as.numeric(dfstrains[i, 3:ncol(dfstrains)])
    prob1 = strainweight[1,]
   
    prob2 = strainweight[2,]
  
    table1 = table(sample(a, size = depth1,replace = TRUE , prob = prob1))
    if (length(table1)==2){
      sequenced$ref1[i] =  table1[["0"]]
      sequenced$alt1[i] =  table1[["1"]]
      
    } else if (names(table1) == "0") {
      sequenced$ref1[i] =  table1[["0"]]
      sequenced$alt1[i] =  0
    } else {
      sequenced$ref1[i] =  0
      sequenced$alt1[i] =  table1[["1"]]
    }
    
    
    table2 =  table(sample(a, size = depth2,replace = TRUE , prob = strainweight[2,]))
    if (length(table2)==2){
      sequenced$ref2[i] =  table2[["0"]]
      sequenced$alt2[i] =  table2[["1"]]
    } else if (names(table2) == "0") {
      sequenced$ref2[i] =  table2[["0"]]
      sequenced$alt2[i] =  0
    } else {
      sequenced$ref2[i] =  0
      sequenced$alt2[i] =  table2[["1"]]
    }
    
  
  }
  sequenced$SNPindL = sequenced$alt1 / sequenced$depth1
  sequenced$SNPindH = sequenced$alt2 / sequenced$depth2
  sequenced$deltaSNP = sequenced$SNPindH - sequenced$SNPindL
  sequenced<-filter(sequenced, depth2>10 &depth1>10)
  
  return(sequenced)
}



############################END OF FUNCTION##############################
#############################START OF CODE###############################
# 
# #35 is the depth
# #SNP value is first probability
# table(rmultinom(35, 1, c(0.54, 0.46) )[1,])
# #number of zeros
# table(rmultinom(35, 1, c(0.54, 0.46) )[1,])[[1]]
# #number of ones
# table(rmultinom(35, 1, c(0.54, 0.46) )[1,])[[2]]
# 

#rmultinom(1, 35, c(0.54, 0.36, 0.1))

# 
# a<- as.numeric(dfstrains2[3, 3:ncol(dfstrains2)])
# table(sample(a, size = depth,replace = TRUE , prob = strainweights2[1,]))

sequenced = sequence(75, SNPindex2, dfstrains2, strainweights2)



