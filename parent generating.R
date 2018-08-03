#hapchr = list of haplotype sizes by chromosomes
#parent1 = first parent (usually KN99) just chr, pos and genotype
#parent2 = second parent
#numstrains = how many strains to generate
genstrain2 <- function(hapchr, parent1, parent2, numstrains = 1) {
  #only include spots where they are different, no NA (no reads), and no alternative 2
  
  compare = data.frame(CHR = parent1$CHR , POS = parent1$POS, p1 = parent1$Genotype, p2 = parent2$Genotype, stringsAsFactors=FALSE)
  compare = filter(compare, p1!=p2 & !is.na(p1) & !is.na(p2) &((p1 == 0)|(p1 ==1)) &((p2 == 0)|(p2 ==1)))
  #makes dataframe with positions of SNPs
  storev2 = compare[1:2]
  storev2$CHR <- factor(storev2$CHR, levels = unique(storev2$CHR))
  for (j in 1:numstrains) {
    # adds empty genotype column for new strain
    storev2[,j+2] = NA
    #tracker of position in dataframe for while loop
    i = 1
    #randomly generates if it should start with parent  one or two
    isParent1= sample (c(TRUE,FALSE), size = 1, replace = TRUE, prob = c(0.5, 0.5))
    #generates if strain should be parental    
    parental = sample(c(TRUE, FALSE), 1,replace=TRUE, prob=c(0.15, 0.85))
    print(parental)
    #if its parental just assign the strain values from only one parent
    if (parental) {
      if(isParent1) {
        storev2[,j+2] = compare[,3]
      }
      else{
        storev2[,j+2] = compare[,4]
      }
    }
    else {
      #for the number of positions
      while (i <= length(storev2$POS)) {
        #generate the length of the haplotype
        length = sample(hapchr[[storev2$CHR[i]]], size = 1, replace = TRUE)
        k = i
        #assign genotype based on parents
        if (isParent1) {
          storev2[i,j+2] = compare[i,3]
       
        } else {
          storev2[i,j+2] = compare[i,4]
          
        }
        i = i+1
        #while the length is less than the target length, and there are still positions left
        while ((storev2$POS[i]-storev2$POS[k]+1)<length & i<=length(storev2$POS)){
          #if the chromosomes are different re sample which parent to use and move on to next segment
          if (storev2$CHR[i]!= storev2$CHR[i-1]){
            isParent1= sample (c(TRUE,FALSE), size = 1, replace = TRUE, prob = c(0.5, 0.5))
            break()
          }
          #otherwise assign the genotype for the correct parent
          else {
            if (isParent1) {
              storev2[i,j+2] = compare[i,3]
            } else {
              storev2[i,j+2] = compare[i,4]
            }
            
            
            
            i = i+1
            
          }
          
        }
        
        
        #switches genotype for the next section
        isParent1 = !isParent1
      }
      
      
      
    }
  }
  
  
  
  #return strains
  return(storev2)
}

#numstrains = how many strains were generated
calcweights2 <- function(numstrains =1) {
  #calculates equal weights for all strains for the starting weights
  strainweights = data.frame(matrix(NA, nrow = 1, ncol = numstrains))
  weight = 1/numstrains
  strainweights[1,] = weight
  return(strainweights)
}



####################################end of functions#########################
####################################start of code/calls######################

#generate parent 1
TDY1993 = allSamples[[2]]

TDY1993 = filter(TDY1993, TDY1993$CHR != "chrM")
#isolates chromosome, position and genotype columns
TDY1993 = TDY1993 %>% dplyr::select(CHR, POS, Genotype)


#generate parent 2
KN99 = allSamples[[1]]

#only includes snps that are filtered

KN99 = filter(KN99, KN99$CHR != "chrM")

#isolates chromosome, position and genotype columns
KN99 = KN99 %>% dplyr::select(CHR, POS, Genotype)




dfstrains2 <- genstrain2(hapchr, KN99, TDY1993, numstrains =   40)

strainweights2 <- calcweights2(numstrains = 40)


