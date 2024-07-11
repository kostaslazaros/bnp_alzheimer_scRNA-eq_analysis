library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(openxlsx)

df <- read.csv("./data/8.data4degs/alzheimer_data_degs_v1.csv")
df$neg_log10_pval <- -log10(df$pvals)
df

# Biostatsquid theme
theme_set(theme_classic(base_size=20) + 
            theme(
              axis.title.y=element_text(face="bold", margin=margin(0,20,0,0), size=rel(1.1), color='black'),
              axis.title.x=element_text(hjust=0.6, face="bold", margin=margin(20,0,0,0), size=rel(1.1), color='black'),
              plot.title=element_text(hjust=0.6)
            ))

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"
df$diffexpressed[df$logfoldchanges > 1 & df$pvals < 0.001] <- "UP"
df$diffexpressed[df$logfoldchanges < -1 & df$pvals < 0.001] <- "DOWN"
head(df[order(df$pvals) & df$diffexpressed == 'DOWN', ])

# For Downregulated genes: prioritize by -log10(pval) first (high significance), then by logfoldchanges (most negative)
top_downregulated <- df[df$diffexpressed == "DOWN", ][order(-df$neg_log10_pval[df$diffexpressed == "DOWN"], df$logfoldchanges[df$diffexpressed == "DOWN"]), ][1:20,]


# For Upregulated genes: prioritize by -log10(pval) first (high significance), then by logfoldchanges (most positive)
top_upregulated <- df[df$diffexpressed == "UP", ][order(-df$neg_log10_pval[df$diffexpressed == "UP"], -df$logfoldchanges[df$diffexpressed == "UP"]), ][1:81,]



# Combine into a list
top_genes_combined <- c(top_downregulated$genes, top_upregulated$genes)
top_genes_combined


df_annotated <- df[df$genes %in% top_genes_combined,]

write.csv(df_annotated, file = "./top100degs.csv", row.names = FALSE)
write.xlsx(df_annotated, file = "./top100degs.xlsx")

p <- ggplot(data = df, aes(x=logfoldchanges, y=neg_log10_pval, col=diffexpressed)) + 
       geom_vline(xintercept=c(-1, 1), col="gray", linetype='dashed') + 
       geom_hline(yintercept= -log10(0.05), col="gray", linetype='dashed') + 
       geom_point(size = 2) + 
       scale_color_manual(values=c("#00AFBB", "grey", "#bb0c00"), labels=c("Downregulated", "Not significant", "Upregulated")) + 
       coord_cartesian(ylim=c(0, 200), xlim=c(-5, 5)) + 
       labs(color='', x=expression("log"[2]*"FC"), y=expression("-log"[10]*"p-value")) + 
       scale_x_continuous(breaks=seq(-10, 10, 2)) + 
       ggtitle("Volcano of DEGs (Disease vs Control)")
  
# Adding labels with ggrepel for better visibility and avoiding overlaps
p + geom_label_repel(
      data=df_annotated,
      aes(label=genes, x=logfoldchanges, y=neg_log10_pval),
      box.padding=0.35,
      point.padding=0.5, 
      size=3, 
      segment.color='grey50', 
      show.legend=FALSE
    )
ggsave("./figures/deg_volcano_v1.png", dpi=600)
  