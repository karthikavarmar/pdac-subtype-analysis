#--------------------------------------------SAM Analysis---------------------------------------


install.packages("siggenes")
library(siggenes)
library(multtest)
install.packages("pamr")
library(pamr)

#expr_rm <- apply(expr_rm, 1, median)

ccpdac.cl <- pca_data[,3]

sam.out <- sam(expr_rm, ccpdac.cl, rand=123)
s <- summary(sam.out) 
capture.output(s, file = "sam_summary.txt")

#plotting the silhouette plot
print(paste0("Plotting Silhouette Plot and Principal Component Analysis Biplot (with #batches and clustering information) to the file ", plotFile))
plotFile = paste0(date_today, "_samplot.pdf")
pdf (plotFile)
plot(sam.out)
dev.off()

#plotting SAM plots for delta = 0.1 (& 32.8) 
plotFile = paste0(date_today, "_samplot_point1.pdf")
pdf (plotFile)
plot(sam.out,0.1)
dev.off()

sam.sum <- summary(sam.out,32.8)
capture.output(sam.sum, file = "samsummary_32point8.txt")

#finding the significant genes
sam.sum@row.sig.genes
sigf<- sam.sum2@mat.sig
capture.output(sigf, file = "siggenes_sam.txt")

#-------------------------------PAM analysis-------------------------------------------------
  
m0 = match(names(summary(sam.out,32.8,entrez=FALSE)@row.sig.genes),rownames(expr_rm))
w0 = which(!is.na(m0)) 
tD = expr_rm[m0[w0],]
#dim(tD)

Ypcen<-tD 
set.seed(123)

data = list(x=as.matrix(Ypcen), y=ccpdac.cl, geneid = rownames(Ypcen), genenames = rownames(Ypcen))
train = pamr.train(data)
cv = pamr.cv(train,data, nfold=5)

pamr.confusion(cv, threshold = cv$threshold[which(cv$error==min(cv$error))][1])

pamcen <- pamr.listgenes(train, data, threshold = cv$threshold[which(cv$error==min(cv$error))][1])

colnames(pamcen)<-c("genes",colnames(pamcen)[-1]); #dim(pamcen)
colnames(pamcen)<-c("genes",colnames(pamcen)[-1]); #dim(pamcen)

pdf(paste0(Sys.Date(),"_ccpdac_sam-pam_32point8_centroids.pdf"))
pamr.plotcv(cv)
dev.off()

write.table(pamcen,paste0(Sys.Date(),"_ccpdac_sam-pam_32point8_centroids.txt"), quote = TRUE, sep = "\t", row.names = FALSE)

capture.output(pamcen, file = "pamcentroids.txt")
