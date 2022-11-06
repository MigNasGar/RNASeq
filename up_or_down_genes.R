#With this code, it's possible to filer UP or DOWN regulated genes from a DESeq2 table
#Author: Miguel A. N. Garcia 
#Laboratory of Molecular Biology of Pathogens (Unifesp/EPM)

install.packages("dplyr")
install.packages("ggrepel")
install.packages("writexl")

library(dplyr)
library(ggplot2)
library(writexl)

#Import the .tab file and defines that there are no header (with header=TRUE files there is some bugs in the analysis)
results <- read.delim('/file.tab', header = FALSE) #insert the pathway where your file is located

#Assigns wich column represents witch data, in the majority of cases the log2FC is in the 2nd column, geneID in the 1st and the adjusted pvalue in the 6st
results <- results %>% mutate(fdr = .[[7]],
                              pvalue = .[[6]],
                              logfc = .[[3]],
                              labels = .[[1]])

#Mekes names for echa legend, down regulated, up regulated or non significant for those genes who do not pass the logical test
down <- unlist(strsplit('Down,Not Sig,Up', split = ","))[1]
notsig <- unlist(strsplit('Down,Not Sig,Up', split = ","))[2]
up <- unlist(strsplit('Down,Not Sig,Up', split = ","))[3]

#Do the logical test for the UP REGULATED genes, with the threshold for the adjusted pvaule = 0.05 
results <- mutate(results, sig = case_when(
                                fdr < 0.05 & logfc > 0.0 ~ up,
                                TRUE ~ notsig))

#Do the logical test for the DOWN REGULATED genes, with the threshold for the adjusted pvaule = 0.05 
results <- mutate(results, sig = case_when(
                                fdr < 0.05 & logfc < -0.0 ~ down,
                                TRUE ~ notsig))

#Shows the table for confirms if there is the up or down genes, with the nonsig 
head(results)

#Writes the table in a Excel file for further analysis or filter spcific genes
write_xlsx(results,"/content\\file_results.xlsx")