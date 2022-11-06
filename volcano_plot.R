#install.packages("dplyr")
#install.packages("ggrepel")
#install.packages("EnhancedVolcano")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

#the .tab archive should not have a header for this to work
results <- read.delim('/file.tab', header = FALSE)

#determinates which data is in wich column, in this case, the data is in the 1,3,6 and 7 column 
results <- results %>% mutate(fdr = .[[7]],
                              pvalue = .[[6]],
                              logfc = .[[3]],
                              labels = .[[1]])

#names for legend
down <- unlist(strsplit('Down,Not Sig,Up', split = ","))[1]
notsig <- unlist(strsplit('Down,Not Sig,Up', split = ","))[2]
up <- unlist(strsplit('Down,Not Sig,Up', split = ","))[3]

#set colours
colours <- setNames(c("cornflowerblue", "grey", "firebrick"), c(down, notsig, up))

#create the significant (sig) column
results <- mutate(results, sig = case_when(
  fdr < 0.05 & logfc > 0.0 ~ up,
  fdr < 0.05 & logfc < -0.0 ~ down,
  TRUE ~ notsig))

#Get top genes by P value
top <- slice_min(results, order_by = pvalue, n = 10)
# Extract into vector
toplabels <- pull(top, labels)
# Label just the top genes in results table
results <- mutate(results, labels = ifelse(labels %in% toplabels, labels, ""))

#saves the plot in a pdf file
pdf("volcano_plot.pdf")
p <- ggplot(data = results, aes(x = logfc, y = -log10(pvalue))) +
  geom_point(aes(colour = sig)) +
  scale_color_manual(values = colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank())

#add labels to the genes
p <- p + geom_text_repel(data = filter(results, labels != ""), aes(label = labels),
                         min.segment.length = 0,
                         max.overlaps = Inf,
                         show.legend = FALSE)

#plot the volcano
p <- p + theme(legend.title = element_blank())
print(p)
dev.off()