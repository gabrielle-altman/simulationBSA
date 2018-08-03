# simulationBSA
BSA simulation code


BSA SIMULATION GUIDELINES:

packages used: ggplot, dplyr, openxlsx


1.	Run program to convert vcf file into allSamples – daniel’s code stored in file haplotypes1strain.R
2.	Upload excel file haplotype_data_chr.xlsx to have haplotype sizes by chromosome (DOWNLOAD openxlsx PACKAGE FOR THIS STEP)
        hapchr <- read.xlsx("haplotype_data_chr.xlsx", sheet = 1)
        for (i in 2:14) {  
          hapchr <- c(hapchr, read.xlsx("haplotype_data_chr.xlsx", sheet = i))
        }
3.	Run “parent generating.R”
4.	Run “selection.R” or “selection2.R
5.	Run “sequencing.R”
6.	Run “bsaPackage.R”




Run “calculating SNP index.R” to calculate SNP index without smoothing or windows



