#dfstrains: the structure holding the strains
#strainweights: beginning and ending weights
calcIndex <- function(dfstrains, strainweights) {
  #create dataframe
  SNPindex <- data.frame(CHR = dfstrains$CHR, POS = dfstrains$POS, indexL = NA, indexH = NA, delta = NA)
  # for (i in 1:nrow(SNPindex)){
  #   sum =0
  #   
  #   for (j in 1:ncol(strainweights)) {
  #     sum = sum+ dfstrains[i,j+2]*strainweights[1,j]
  #   }
  #   SNPindex$indexL[i] = sum
  # }
  
  
  for (i in 1:nrow(SNPindex)) {
    #calculate low bulk
    SNPindex$indexL[i] = sum(strainweights[1,]*dfstrains[i,3:ncol(dfstrains)])
    #calculate high bulk
    SNPindex$indexH[i] = sum(strainweights[2,]*dfstrains[i,3:ncol(dfstrains)])
    #calculate deltaSNP
    SNPindex$delta[i]= SNPindex$indexH[i] - SNPindex$indexL[i]
  }
  
  return(SNPindex)
}

####################################end of function##################################

SNPindex2 <- calcIndex(dfstrains2, strainweights2)
#hist(SNPindex2$indexL)


#ggplot(data=SNPindex2)+geom_line(aes(x=POS,y=delta))+facet_wrap(~CHR)


#ggplot(data=copy)+geom_point(aes(x=POS,y=V43, color = V3), size = 0.5) + geom_point(aes(x=POS,y=V44, color = V4), size = 0.5)+geom_point(aes(x=POS,y=V45, color = V5), size = 0.5)+geom_point(aes(x=POS,y=V46, color = V6), size = 0.5)+facet_wrap(~CHR)



