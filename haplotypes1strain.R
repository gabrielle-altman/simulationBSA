library(xlsx)
library(ggplot2)
library(stringr)
library(knitr)
library(ggrepel)
library(dplyr)
library(reshape2)

#dplyr and subset function

setwd("~/Documents/lab work")

#This one is for the TDY1993 analysis.
tab=read.table("TDY1993.R.vcf", as.is=T, header=T)

##########################################################################################################################################

tabler=function(i=1, minDepth=5, cutoff=0.75){
  c=(str_split_fixed(tab[,9+i],":",n=8)) #This splits the FORMAT itens from the vcf into columns.
  colnames(c)=c("Genotype_GT","Depth_DP","AD","Reference_RO","QR","Alternative_AO","QA","GL")
  
  d=str_split_fixed(c[,3], ",",n=3) #This splits the two alternative allele counts into different columns.
  d=apply(d,2,as.numeric)
  e=str_split_fixed(tab[,5], ",",n=3)
  for(j in 1:nrow(d)){ #This for loop turns all NAs in the both columns into 0 or ".".
    if(is.na(d[j,3])){d[j,3]=0}
    if(is.na(d[j,2])){d[j,2]=0}
    if(is.na(d[j,1])){d[j,1]=0}
    if(e[j,3]==""){e[j,3]="."}
    if(e[j,2]==""){e[j,2]="."}
    if(e[j,1]==""){e[j,1]="."}
  }
  
  data2=data.frame(CHR=tab[,1], POS=tab[,2], Ref_allele= tab[,4], ALT1=e[,1],ALT2=e[,2],ALT3=e[,3], Genotype=as.numeric(c[,1]),
                   Depth=as.numeric(c[,2]),Reference= as.numeric(c[,4]), Alternative1=d[,2], Alternative2=d[,3], 
                   filtGenotype=rep(".",nrow(e)), stringsAsFactors = FALSE) 
  data2$filtGenotype=as.character(data2$filtGenotype)
  for(j in 1:nrow(data2)){ ## Loop for defining the filtered Genotype. Filter by minimum depth and Genotype cutoff. The cutoff is the
    # minimum percentage of reads necessary to define a genotype. Default = 0.75.
    if(is.na(data2$Reference[j])){data2[j,12]="noReads"}
    else if(data2[j,8]< (minDepth+1)){data2[j,12]="lowDepth"}
    else if(data2[j,8]> minDepth & data2$Reference[j]>cutoff*data2$Depth[j]){data2[j,12]="Reference"}
    else if(data2[j,8]>minDepth & data2$Alternative1[j]>cutoff*data2$Depth[j]){data2[j,12]="Alternative1"}
    else if(data2[j,8]>minDepth & data2$Alternative2[j]>cutoff*data2$Depth[j]){data2[j,12]="Alternative2"}
    else if(data2[j,8]>minDepth & (data2$Alternative2[j]<=cutoff*data2$Depth[j] &  data2$Alternative1[j]<=cutoff*data2$Depth[j] &
                                   data2$Reference[j]<=cutoff*data2$Depth[j])  ){data2[j,12]="undefinedGenotype"}
    
  }
  
  data2$Ref_percentage=data2$Reference/data2$Depth
  data2$Alt1_percentage=data2$Alternative1/data2$Depth
  data2$Alt2_percentage=data2$Alternative2/data2$Depth
  rm(c)
  rm(d)
  return(data2)
}


##########################################################################################################################################
##CODE STARTS HERE
tab1=tab[tabler()$Genotype==0,]

allSamples=list()

# These two objects will be used for calling dataframes in the list.
name_cols=colnames(tab[,10:(ncol(tab))])

Strain=name_cols[order(name_cols)] # These two objects will be used for calling dataframes in the list.


#This will fill the allSamples list object with the table (from tabler function) from each strain as each element. 
for (z in 1:length(name_cols)){ 
  allSamples[[name_cols[z]]] <- tabler(i=z)
  
}
#########################RUN UNTIL HERE FOR STEP ONE#####################
#allSamplesc = allSamples
#allSamplesc = filter(allSamplesc, allSamplesc[[4]][7]!=1)

#isolates just this strain
TDY3190 = allSamples[[4]]

#only includes snps that are filtered
TDY3190 = filter(TDY3190, filtGenotype == "Reference" | filtGenotype == "Alternative1")
#isolates chromosome, position and genotype columns
TDY3190 = TDY3190 %>% dplyr::select(CHR, POS, Genotype)

#initializes list of breakpoints
brkpts <- list("chr1" =1, "chr2" =1, "chr3" =1, "chr4" =1, "chr5" =1, "chr6" =1, "chr7" =1, "chr8" =1, "chr9" =1, "chr10" =1, "chr11" =1, "chr12" =1,"chr13" = 1, "chr14" = 1)
#list of chromosome lengths
chrlen = c(2291500, 1621676,1574972,1084805,1814975,1422463,1399209,1398693,1186813,1059962,1562107,774060,756017,942472)

#filters out anything that is NA
TDY3190= filter(TDY3190, !is.na(Genotype ))


#for one less than the number of SNPs
for (i in 1:(nrow(TDY3190)-1)) {
  #if the genotype and the next one are the same move on
 if (TDY3190[i,3]==TDY3190[i+1,3]){
    next()
 }
  #if it is on chrM skip it
  else if(TDY3190[i,1]=="chrM"){
    next()
  }
  #if the SNPs are on the same chromosome
  else if (TDY3190[i,1]==TDY3190[i+1,1]){
    #set the brkpt equal to the average of the two positions
    brkpts[[TDY3190[i,1]]]=c(brkpts[[TDY3190[i,1]]], floor((TDY3190[i,2]+TDY3190[i+1,2])/2))
  }
  #otherwise move on
  else {
    next()
  }
}

#add in the end of the chromosome
for (j in 1:length(chrlen)) {
  brkpts[[j]] = c(brkpts[[j]], chrlen[[j]])
}


#initialize haplotype dataframe
haps <- data.frame(chr1= numeric() , chr2= numeric() , chr3= numeric() , chr4= numeric(), chr5= numeric() , chr6= numeric(), chr7= numeric() , chr8= numeric(), chr9= numeric() , chr10= numeric(), chr11= numeric() , chr12= numeric(), chr13= numeric() , chr14= numeric())

#for each chromosome
for (j in 1:ncol(haps)) {
  #for the number of breakpoints minus one
  for (i in 1:(length(brkpts[[j]])-1)) {
    #calculate length btw breakpoints and add to dataframe
    haps[i,j]=  brkpts[[j]][i+1]-brkpts[[j]][i]+1
    
  }
  
}

#list all lengths together
haplist = c (haps[[1]], haps[[2]], haps[[3]], haps[[4]], haps[[5]], haps[[6]], haps[[7]], haps[[8]], haps[[9]], haps[[10]], haps[[11]], haps[[12]], haps[[13]], haps[[14]])
#remove NA
haplist = haplist[!is.na(haplist)]
#histogram to visualize cutoff sizes
hist(haplist, breaks = 20)
histinfo = hist(haplist, breaks = 40)
