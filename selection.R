
#m = number of mutations to select
#strains = "Sequences" of strains
#weights = proportions of each strain
#compare: filtered parents
genmutation<- function(m, strains, weights, mean, sd, compare){
  #randomly select m mutations
  mutPOS = sample(1:length(strains$POS), size = m ) 
  #make dataframe with mutations
  mutation= strains[mutPOS,]
  #add column for "selection coefficients"
  mutation[, (length(mutation)+1)]=NA
  names(mutation)[length(mutation)] = "coef"
  #sample coefficients from normal distrubution
  for (i in 1:length(mutation$coef)) {
    mutation$coef[i]=sample(rnorm(length(compare[[1]]), mean  = mean, sd = sd) , size = 1)
  }
  #set length of the row that holds the coefficient of each strain
  lenr = length(mutation[[1]])+1
  
  for (i in 3:(length(mutation)-1)) {
    #multiply genotypes in one strain by their coefficients and sum for total score of the strain
    mutation[ lenr, i] = sum(mutation[1:(lenr-1),i]*mutation[1:(lenr-1),length(mutation)])
    #convert from log2 fold change 
    mutation[ lenr, i] = 2^(mutation[ lenr, i])
  }
  #adjust the weights
  weights[2,] = weights[1,]*mutation[ lenr, 3:(length(mutation)-1) ]
  return(list(mutation,weights))
}



############################END OF FUNCTION##############################
#############################START OF CODE###############################


results = genmutation(50, dfstrains2, strainweights2, -0.1, 1, compare)
mutations = results[[1]]
strainweights2 = results[[2]]
#normalize to 1
strainweights2[2,] = strainweights2[2,]/(sum(strainweights2[2, ]))











