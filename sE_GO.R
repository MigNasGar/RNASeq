#simplifyEnrichment code for plot heat map graphs using a list of GO IDs

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("simplifyEnrichment") #package installation is only required if you have never used it before

#######

library(simplifyEnrichment)
set.seed(888)
go_id = (Column1 = c("GO:0006139","GO:1901360","GO:0006457","GO:0006725","GO:0046483","GO:0019693",
                     "GO:0009259",	"GO:0071840",	"GO:0009117",	"GO:0006753",	"GO:0009150",	"GO:1902600",
                     "GO:0044281",	"GO:0046034",	"GO:0008150",	"GO:0055086",	"GO:0006163",	"GO:0072521",
                     "GO:0007088",	"GO:0007052",	"GO:0042026",	"GO:0019725",	"GO:0043094",	"GO:0006873",	
                     "GO:0016053")) #data frame of GO IDs as example, it is ideal to use at least >20 IDs for a good plot and analysis							

mat = GO_similarity(go_id) #uses the provided dataframe for the analysis

GO_similarity(go_id, measure = "Wang") #Here I use de "Wang" method for the analysis, bur other methods can be use in accordance with the package documentation

df = simplifyGO(mat) #plot the heatmap
