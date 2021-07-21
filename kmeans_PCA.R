library("readxl") #Reads excel-files
library(dplyr) #Data manipulation
library(tibble) #Data frame manipulation
library(ggplot2) #Data visualization
library(readr) #CSV file I/O, e.g. the read_csv function

#Reading PDAC samples 
print("Reading GExp PDAC data :")
ccpdac <- read.delim('C:/Users/Karthika/Documents/Anguraj/01Clustering/Sample/2019-04-05_2019-04-05_92_COMBAT_corrected_rnaseq_PDAC_4_batches_counts_TMM_CPM_log2_ENS_ID_removed_80_55968_with_gene_name_data_sd0.txt', stringsAsFactors = FALSE)
#grepl(" /// ", ccpdac$Gene)

ccpdac.genes<-data.frame(Gene='', row='', check.names = FALSE)

#splitting gene symbols in certain rows
split_gene_symbols<- function(row.num){
  
  #getting the gene symbol associated with a particular row
  string <- ccpdac$Gene[row.num]
  
  #splitting the genes (if there are multiple gene names)
  genes <- unlist( strsplit(string, split=" /// ", fixed=TRUE))
  
  #add all genes to a list
  lapply((1:length(genes)), function(x){
    c(gene=genes[x], row=row.num)
    })
}

#Adds all genes from all rows to a list
ccpdac.genes<-t(data.frame(lapply(c(1:length(ccpdac$Gene)), split_gene_symbols)))

row.names(ccpdac.genes)= c(1:dim(ccpdac.genes)[1])

#Reading the Gene list-annotated data
geneset <-read_excel("./20210420_144_Gene_PDACMarvel_KD.xlsx")

#Matching the genes from the dataset with the geneset
genes<-as.character(as.matrix(geneset)[,1])
id <- as.character(ccpdac.genes[,1])
mat<- match (genes, id)
dim(mat)

colsum <- (colSums(ccpdac[,-1]))

expr<- ccpdac[na.omit(ccpdac.genes[mat, 2]), ]
write.csv(expr, "exprwithgenes.csv", row.names= TRUE)

#removing the first column containing gene symbols before clustering
expr_rm <- expr[,-1]
write.csv(expr_rm, "exprdat.csv", row.names= TRUE)


#-----------------------------------K_Means Clustering-----------------------------------


expr_dat <- t(expr_rm)
write.csv(expr_dat, "exprdat.csv", row.names= TRUE)

setwd("C:/Users/Karthika/Documents/Anguraj/01Clustering/Sample/PDAC_actual")

#calculating distance matrix
distMatrix <- dist(expr_dat, method = "euclidean")

#using the silhouette method for determining the optimal number of clusters for k-means

#clustering data for various k
cluster_data<- lapply(c(2:7), function(x){kmeans(expr_dat,
                                                 centers= x,
                                                 iter.max = 1000,
                                                 nstart=25,
                                                 set.seed(12345))})


silh <- as.data.frame(lapply(c(2:7),  function(x){
  c(x, mean(cluster::silhouette(cluster_data[[x-1]]$cluster, distMatrix)[,3]))}))

colnames(silh)<- c(2:7)
row.names(silh) <- c("k", "silWidth")

silh <- as.data.frame(t(silh))

#writing avg sil width to file
date_today <- as.character(format(Sys.Date(), "%Y%m%d"))
write.table(silh, file = paste0(date_today, "_ccpdac_filtered_avg_silhouette_width_k-means.txt"),
            sep = "\t")


#finding the k with maximum avg. silhouette width
max.sil <- silh[1+which(silh[2:dim(silh)[1],2]==max(silh[2:dim(silh)[1],2])),]
print(paste0("k=", max.sil[1], " is optimal with Average Silhouette Width = ", max.sil[2]))
opt.k <- as.numeric(max.sil[1])

#writing the clustering information to file
clustered_data<-cluster_data[[opt.k - 1]]
cluster.info <- cbind(names(clustered_data$cluster), clustered_data$cluster)
colnames(cluster.info)<- c("Samples", "Cluster")
write.table(cluster.info,paste0( date_today, "_ccpdac_filtered_k-means_cluster_info.txt"),
            sep = "\t", row.names = FALSE)


#calculating Principal Components
pca.data <- prcomp(expr_dat, center = TRUE, scale. = TRUE)

#matching sample IDs from PCA and cluster info
id<- row.names(pca.data$x)
sample.match <- match(id, row.names(cluster.info))

#combining PCs, class information, K-means cluster information
pca_data <- data.frame(pca.data$x[, 1:2], cluster.info[sample.match, 2])
colnames(pca_data)[3] <- "cluster"

pca_data[,3] <- as.factor(pca_data[,3])

#plotting silhouette plot and PCA plot with batches and kmeans information for opt.k
plotFile = paste0(date_today, "_ccpdac_filtered_pca_k-means.pdf")

print(paste0("Plotting Silhouette Plot and Principal Component Analysis biplot (with batches and clustering information) to the file ", plotFile))
pdf (plotFile)



#silhouette plot
silh.plot <- ggplot2::ggplot(data = silh, ggplot2::aes(x =k, y =silWidth)) +
  ggplot2::geom_vline(ggplot2::aes(xintercept = opt.k), size =1,
                      linetype = 'dashed', colour = "dodgerblue3")+
  ggplot2::geom_point(size = 2.5, colour = "orange")+
  ggplot2::geom_path(size = 1.25, colour = "orange") +
  ggplot2::labs(x = "Number of clusters, k",
                y = "Average Silhouette Width",
                title = "Optimal number of clusters for k-means") +
  
  ggplot2::theme_classic()

plot(silh.plot)



#PCA (with batch) and kmeans info plot
kmeans_pca_plot <- ggplot2::ggplot(data = pca_data) +
  
    ggplot2::geom_point(ggplot2::aes(x=PC1,
                                      y=PC2,
                                      colour =cluster),
                                      size=6, alpha = 0.6)+
    ggplot2::labs(title= paste0("PCA plot for k = ", opt.k),
                colour = "k-means cluster")+
    ggplot2::theme(
    
      #Adjusting axis titles, lines and text
      axis.title = ggplot2::element_text(size = 15),
      axis.line = ggplot2::element_line(size =0.75, colour = "black"),
      axis.text = ggplot2::element_text(size=15, colour ="black"),
    
      #Center align the title
      plot.title = ggplot2::element_text(face = "bold", hjust =0.5, size =20),
    
      #Adjust legend title and text
      legend.title = ggplot2::element_text(size = 15, face = "bold"),
      legend.text = ggplot2::element_text(size = 15),
    
      # Remove panel border
      panel.border = ggplot2::element_rect(fill=NA, size= 0.75),
    
      # Remove panel grid lines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
    
      # Remove panel background
      panel.background = ggplot2::element_blank()) +
  
    ggplot2::scale_color_manual(values=c("#fcba03",  "#19e01c", "#ff470f",
                                       "#0fdbca", "#ff217e","#405ce6",
                                       "#6b6769","#b264ed"))

plot(kmeans_pca_plot)

dev.off()
