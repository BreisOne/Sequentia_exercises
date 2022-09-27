libraries <- c('DESeq2','ggVennDiagram', 'cowplot','reactome.db','org.Hs.eg.db','ggpubr','rstatix','readxl','dendsort','apeglm','VennDiagram', 'RColorBrewer', 'pheatmap','tidyverse','ggrepel')
lapply(libraries,library, character.only = TRUE)

setwd("C:/Users/Brise/OneDrive - Universidade de Vigo/Sequentia/E3")

#####Functions

save_plot <- function(plot,name, w, h){
  pdf(file=paste0("./Figuras",name), width=w, height=h)
  plot
}

volcano_plot <- function(df, title, colors){
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = interaction(threshold_padj,Overexpressed,Underexpressed), label = row.names(df))) + 
    geom_point() +
    scale_colour_manual(values = colors)+
    #geom_label_repel(data = head(df,20), aes(label=GS), max.overlaps = Inf) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#AA4499")+
    geom_vline(xintercept = 1.5, linetype="dashed", color = "#AA4499")+
    geom_vline(xintercept = -1.5, linetype="dashed", color = "#AA4499")+
    scale_x_continuous(limits = c(-3, 3), breaks = seq(from = -3, to = 3, by = 0.5))+
    scale_y_continuous(limits = c(0, 14), breaks = seq(from = 0, to = 14, by = 4))+
    ggtitle(title)+
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") + 
    theme_bw()+
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))
}

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

##Import data

raw_data <- read.delim("counts.txt")

AT_rownames <- raw_data[,1]

raw_data <- raw_data[,c(8:9,12,10:11,7)]

raw_data <- raw_data  + 1

row.names(raw_data) <- AT_rownames

 

# Generar el dataframe del diseño experimental.

l <- colnames(raw_data)
c <- as.factor(c(rep("CTRL",3), rep("STR",3)))
r <- as.factor(c(c(rep(seq(1,3),2))))
colnames <- c("condition","replicate")
exp_design <- data.frame(c,r, stringsAsFactors = FALSE)
colnames(exp_design) <- colnames
rownames(exp_design) <- l


#Create the DESeq2 object
AT_DESeq <- DESeqDataSetFromMatrix(countData = raw_data,
                                   colData = exp_design,
                                   design = ~condition)

####Control de calidad de las muestras####

#Determine the size factors to use for normalization
AT_DESeq <- estimateSizeFactors(AT_DESeq)

#Transform the counts matrix
AT_DESeq_VS <- vst(AT_DESeq, blind=TRUE)
AT_DESeq_rlog <- rlog(AT_DESeq, blind=FALSE)

#Extract the matrix of transformed counts
vsd_mat_AT <- assay(AT_DESeq_VS)
rlog_mat_AT <- assay(AT_DESeq_rlog)
#Compute the correlation values between samples
vsd_cor_AT <- cor(vsd_mat_AT) 

# Plot the heatmap correlation
pheatmap(vsd_cor_AT, annotation = select(exp_design, condition))

# Plot the PCA of PC1 and PC2
AT_PCA <- plotPCA(AT_DESeq_VS, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(AT_PCA, "percentVar"))
ggplot(AT_PCA, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=5) +
  ggtitle("PCA plot - RNA-seq AT genes")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()+
  theme(
    axis.title=element_text(size=12,face="bold")
  )

####Differential expression analysis####
AT_Analysis <- DESeq(AT_DESeq)

AT_res <- results(AT_Analysis, 
                  contrast = c("condition","STR","CTRL"),
                  alpha = 0.05,
                  lfcThreshold = 0)

AT_res_all <- data.frame(AT_res)

AT_res_all <- AT_res_all[order(AT_res_all$padj),]

AT_res_all <- AT_res_all %>% mutate(threshold_padj = padj < 0.05,Overexpressed = log2FoldChange >= 1.5,Underexpressed = log2FoldChange <= -1.5)

AT_res_all <- AT_res_all[complete.cases(AT_res_all),]

AT_res_sig <- AT_res_all[AT_res_all$threshold_padj==TRUE&c(AT_res_all$Overexpressed==TRUE|AT_res_all$Underexpressed==TRUE),]

##Heatmap

heat_colors <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)

colours_heatmap_annotation = list(
  condition = c(CTRL = "#6395ED", STR = "#FF4500"),
  replicates = c("1" = "#FACBB3", "2" = "#B3CDE3", "3" = "#CCEBC5"))

#Create dendrogram for samples
mat_cluster_cols <- hclust(dist(t(rlog_mat_AT[row.names(AT_res_sig),])))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")

#Reverse de order in the dendrogram for samples
sort_hclust <- function(...) as.hclust(rev(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")


heat <- pheatmap(rlog_mat_AT[row.names(AT_res_sig),],
                 color = heat_colors,
                 cluster_rows = TRUE,
                 cluster_cols = mat_cluster_cols,
                 show_rownames = FALSE,
                 legend_breaks = c(-1.9,-1,0,1,1.9),
                 annotation_col = select(exp_design, replicate, condition),
                 annotation_colors = colours_heatmap_annotation,
                 scale = "row",
                 main = "RNA-seq AT STR vs CTRL"
)
heat

##Volcano-plot

AT_volcanoplot <-volcano_plot(AT_res_all,"AT STR vs AT CTRL", c("gray","gray","gray","#DC3220","gray","#78E02C"))

AT_volcanoplot

write.csv2(AT_res_sig, "AT_res_sig_DEG.csv", row.names = TRUE)

####GO terms Overlapping genes
GO_terms <- read_excel("Goterms_Analysis.xlsx", 
                                   sheet = "Go Terms")

GO_terms <- GO_terms[GO_terms$Log10_FDR>12,]
#GO_terms_overlapping$`GO type` <- factor(GO_terms_overlapping$`GO type`, levels = GO_terms_overlapping$`GO type`[order(data3$y, decreasing = TRUE)]))
GO_terms<- GO_terms %>% 
                mutate(ordering = as.numeric(factor(Category, levels = c("CC","MF","BP")))+GO_terms$FDR*10000000000) 

ggplot(
  GO_terms,
  aes(
    x=reorder(Terms, desc(ordering)),
    y=`Log10_FDR`, 
    fill= Category, 
    #group=Term,
    #label=Term
  )
) +
  geom_bar(
    stat="identity",
    color="black",
    position="dodge",
    alpha=.8, 
    width=.9
  ) +
  ggtitle("GO terms for overlapping genes") +
  theme_classic()+
  theme()+
  coord_flip()+
  xlab(" ")+
  geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1)
