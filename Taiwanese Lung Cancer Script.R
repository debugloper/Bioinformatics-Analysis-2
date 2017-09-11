# Use of several multivariate analysis of microarray data analysis

# load bioconductor libraries
source("http://bioconductor.org/biocLite.R")
biocLite("affy"); biocLite();biocLite("GEOquery")
library(affy);library(GEOquery);library(Biobase);
rawdata <- ReadAffy() ## call all CELL files

## Normalization of gene expression ##
expr <- rma(rawdata) ##Robust Multichip Average
rmadata <- data.frame(exprs(expr)) ; colnames(rmadata)<-substr(colnames(rmadata), 1, 9)

##### patient information ######
gds3837 <- getGEO('GDS3837', destdir="."); gds3837<- getGEO(filename='GDS3837.soft.gz')
patientinfo <- Columns(gds3837)[,1:5] ## sample information
patientinfo$age <- as.numeric(substr(patientinfo$age, 1, 2))
patientinfo <- patientinfo[order(patientinfo$sample),]
patientinfo$Tumorstage <- factor(substr(patientinfo$other, 14, 14))

###### Matrix of difference between normal and cancer patients expression ##########
diffdf <- matrix(ncol=60,nrow=nrow(rmadata));
for ( i in 1:60) { for ( j in 61:120) {
	if (patientinfo$individual[i] == patientinfo$individual[j]) {
		x <- as.character(patientinfo$sample[i]); y <- as.character(patientinfo$sample[j]);
		diff <- rmadata[,x] - rmadata[,y];
diffdf[,i] <- as.matrix(diff,ncol=1)
} }
}

colnames(diffdf)[1:60]<-paste(patientinfo[,5][1:60]); rownames(diffdf)<-rownames(rmadata)

# tumor & control paired T test
Ttesting <- function(x) {
	tumor <- x[patientinfo$disease.state == "lung cancer"]
	control <- x[patientinfo$disease.state == "control"]
	t_test <- t.test(tumor, control, paired = TRUE)
	p <- t_test$p.value; foldchange <- t_test$estimate; c(p, foldchange)
}

results <- t(apply(rmadata[-(54614:54675),],1,Ttesting))
rownames(results) <- rownames(rmadata[-(54614:54675),]); colnames(results) <- c("p","logFC")
indexy <- which(-log10(results[,"p"])>16) ## only if p-value < 1.0e-16
indexx1 <- which(-log10(results[,"p"])>16&results[,"logFC"]< 0) ## logFC < 0, p-value< 1.0e-16
indexx2 <- which(-log10(results[,"p"])>16&results[,"logFC"]> 0) ## logFC > 0, p-value< 1.0e-16

########## volcano plot #############
plot(results[,"logFC"],-log10(results[,"p"]),xlab="log2 fold change", ylab="-log10(p-value)",
pch=16,cex=0.6); abline(h=16,lty=3,col=6);
points(results[,"logFC"][indexx1], -log10(results[,"p"])[indexx1],col="green",pch=16,cex=0.6)
points(results[,"logFC"][indexx2], -log10(results[,"p"])[indexx2],col="red",pch=16,cex=0.6)

#### Principal Components Analysis ######
pca<-prcomp(t(rmadata[indexy,]))
plot_col <- (1:length(levels(patientinfo$disease.state)))[as.numeric(patientinfo$disease.state)]
plot(pca,type="line",main="screeplot"); plot(pca$x,col=c(2,4)[plot_col],pch=plot_col+14,cex=.8)
text(pca$x[,1]+0.5,pca$x[,2]+0.5,patientinfo$individual,cex=.5)
for (i in 1:60)segments(pca$x[i,1],pca$x[i,2],pca$x[i+60,1],pca$x[i+60,2],cex=.05,lty=3)
legend(10,-15,c("Normal", "Cancer"),pch=c(15,16),col=c(2,4), bty = "n")

#### Independent Component Analysis ####
library(fastICA)
ica<-fastICA(t(rmadata[indexy,]),2)
plot(ica$S,col=c(2,4)[plot_col],pch=plot_col+14,cex=.8,xlab="ICA comp1",ylab="ICA comp2")
text(ica$S[,1]+0.005,ica$S[,2]+0.05,patientinfo$individual,cex=.5)
for (i in 1:60)segments(ica$S[i,1],ica$S[i,2],ica$S[i+60,1],ica$S[i+60,2],cex=.05,lty=3)
legend(-3,-1,c("Normal", "Cancer"),pch=c(15,16),col=c(2,4))

########### Hierarchical clustering on the whole patients
gene <- diffdf[indexy,]; ### gene with p-value < 1.0e-16

library(gplots)
col1 <- c(seq(-4,-1,length=100),seq(-1,1,length=100),seq(1,4,length=100))
palette <- colorRampPalette(c("green", "black", "red"))
heatmap.2(gene, col=palette, breaks=col1, density.info="none", trace="none",
symm=F,symkey=F,symbreaks=T, scale="none")

########### Hierarchical clustering on 3 tumorstages
patientinfo[which(patientinfo[,"Tumorstage"]==4),"Tumorstage"] <- 3
genediff <- diffdf[-(54614:54675),]
tumor <- factor(patientinfo[,"Tumorstage"][1:60]); p_value <- NULL

# ANOVA for the whole genes
for (i in 1:nrow(genediff)){ p_value[i] <- anova(lm(genediff[i,]~tumor))$Pr[1] }
genediff <- cbind(genediff,p_value)

# Working out average per tumor stage for genes with p-value < 1.0e-16
y <- genediff[indexy,]; temp <- NULL; tummean <- NULL;
for (j in 1:699){ temp <- sapply(split(y[j,],tumor),mean); tummean <- rbind(tummean,temp)}
rownames(tummean) <- rownames(y)

# Drawing a heatmap for genes with p-value < 1.0e-16 by diving by tumor stage
heatmap.2(tummean, col=palette, breaks=col1, density.info="none", trace="none",
symm=F,symkey=F,symbreaks=T, scale="none")

library(randomForest) ##### RandomForest clustering

###### Unsupervised randomForest function ########
#### Steve Horvath and Tao Shi ######
collect.garbage <- function(){ while (gc()[2,4] != gc()[2,4]){} }

RFdist <- function(datRF, mtry1, no.tree, no.rep, addcl1 = T, addcl2 = T, imp = T, oob.prox1 = T, proxConver = F) {
	synthetic1 <- function(dat) {
		sample1 <- function(X) { sample (X, replace = T)}
		g1 <- function(dat) {apply(dat, 2, sample1) }
		nrow1 <- dim(dat)[[1]]; yy <- rep(c(1,2), c(nrow1, nrow1));
		data.frame(cbind(yy, rbind(dat, data.frame(g1(dat)))))
	}

	synthetic2 <- function(dat) {
		sample2 <- function(X) { runif(length(X), min=min(X), max=max(X))}
		g2 <- function(dat) { apply(dat,2,sample2)}
		nrow1 <- dim(dat)[[1]]; yy <-rep(c(1,2), c(nrow1, nrow1));
		data.frame(cbind(yy, rbind(dat, data.frame(g2(dat)))))
	}

	cleandist <- function(x) {
		x1 <- as.dist(x); x1[x1<=0] <- 0.000000001; as.matrix(x1)
	}

	nrow1 <- dim(datRF)[[1]]; ncol1 <- dim(datRF)[[2]]
	RFproxAddcl1 <- matrix(0,nrow=nrow1,ncol=nrow1)
	RFproxAddcl2 <- matrix(0,nrow=nrow1,ncol=nrow1)
	RFprox1Conver <- cbind(1:no.rep,matrix(0,(no.rep),3))
	RFprox2Conver <- cbind(1:no.rep,matrix(0,(no.rep),3))
	RFimportance1 <- matrix(0, nrow=ncol1, ncol=4)
	RFimportance2 <- matrix(0, nrow=ncol1, ncol=4)
	RFerrrate1 <- 0; RFerrrate2 <- 0
	rep1 <- rep(666,2*nrow1)

	#adcl1
	if (addcl1) {
		for (i in c(0:no.rep)) {
			index1 <- sample(c(1:(2*nrow1))); rep1[index1] <- c(1:(2*nrow1))
			datRFsyn <- synthetic1(datRF)[index1,]; yy <- datRFsyn[1,];
			RF1 <- randomForest(factor(yy)~., data =datRFsyn[,-1], ntree = no.tree, oob.prox = oob.prox1,
			proximity = TRUE, do.trace = F, mtry = mtry1, importance = imp)
			collect.garbage(); RF1prox <- RF1$proximity[rep1, rep1]

			if (i > 0) {
				if (i > 1) {
					xx <- ((RFproxAddcl1 + (RF1prox[c(1:nrow1), c(1:nrow1)]))/i) - (RFproxAddcl1/(i-1))
					yy <- mean(c(as.dist((RFproxAddcl1 + (RF1prox[c(1:nrow1), c(1:nrow1)]))/i)))
					RFprox1Conver[i,2] <- max(abs(c(as.dist(xx))))
					RFprox1Conver[i,3] <- mean((c(as.dsit(xx)))^2); RFprox1Conver[i,4] <- yy
				}
				RFproxAddcl1 <- RFproxAddcl1 + (RF1prox[c(1:nrow1), c(1:nrow1)])
				if (imp) { RFimportance1 <- RFimportance1 + 1/no.rep*(RF1$importance)}
				RFerrate1 <- RFerrate1 + 1/no.rep*(RF1$err.rate[no.tree])
			}

		}
	}

	# addcl2

	if (addcl2) {
		for (i in c(0:no.rep)) {
			index1 <- sample(c(1:(2*nrow1))); rep1[index1] <- c(1:(2*nrow1))
			datRFsyn <- synthetic2(datRF)[index1,]; yy <- datRFsyn[,1]
			RF2 <- randomForest(factor(yy)~., data = datRFsyn[,-1], ntree = no.tree, oob.prox = oob.prox1,
				proximity = TRUE, do.trace = F, mtry = mtry1, importance = imp)

			collect.garbage(); RF2prox <- RF2$proximity[rep1, rep1]

			if (i > 0) {
				if (i > 1) {
					xx <- ((RFproxAddcl2 + (RF2prox[c(1:nrow1), c(1:nrow1)]))/i) - (RFproxAddcl2/(i-1))
					yy <- mean(c(as.dist((RFproxAddcl2 + (RF2prox[c(1:nrow1), c(1:nrow1)]))/i)))
					RFprox2Conver[i,2] <- max(abs(c(as.dist(xx))))
					RFprox2Conver[i,3] <- mean((c(as.dist(xx)))^2); RFprox2Conver[i,4] <- yy
				}

				RFproxAddcl2 <- RFproxAddcl2 + (RF2prox[c(1:nrow1), c(1:nrow1)])
				if (imp) { RFimportance2 <- RFimportance2 + 1/no.rep*(Rf2$importance)}
				RFerrate2 <- RFerrate2 + 1/no.rep*(RF2$err.rate[no.tree])
			}
		}
	}

	distRFAddcl1 <- cleandist(sqrt(1-RFproxAddcl1/no.rep))
	distRFAddcl2 <- cleandist(sqrt(1-RFproxAddcl2/no.rep))
	distRF <- list(cl1 = NULL, err1 = NULL, imp1 = NULL, prox1Conver = NULL,
	cl2 = NULL, err2 = NULL, imp2 = NULL, prox2Conver = NULL)

	if (addcl1) {
		distRF$sc1 <- distRFAddcl1; distRF$err1 <- RFerrrate1
		if (imp) distRF$imp1 <- RFimportance1
		if (proxConver) distRF$prox1Conver <- RFprox1Conver
	}

	if (addcl2) {
		distRF$cl2 <- distRFAddcl2; distRF$err2 <- RFerrate2
		if (imp) distRF$imp2 <- RFimportance2
		if (proxConver) distRF$prox2Conver <- RFprox2Conver
	}

	distRF
}

## Unsupervised randomForest ( RandomForest clustering ) ##
dataRF <- as.data.frame(gene) ; names(dataRF) <- paste("Sample",names(dataRF),sep="")
clusteringRF <- RFdist(dataRF, mtry1 = 3, 1000, 100, addcl1 = T,addcl2 = T,imp = T, oob.prox1 = T)
cmd1 <- cmdscale(as.dist(clusteringRF$cl1),2); cmd2 <- cmdscale(as.dist(clusteringRF$cl2),2)

g1 <- NULL; g1[which(-log10(results[,"p"]) > 16&results[,"logFC"] < 0)] < -1 ## down
g1[which(-log10(results[,"p"]) > 16&results[,"logFC"] > 0)] <- 2 ##up
g1 <- na.omit(g1); g1 <- as.factor(g1);

plot(cmd1,xlab="Scaling Dimension1",ylab="Scaling Dimension2",col = c(3,2)[g1],pch = 16)
plot(cmd2,xlab="Scaling Dimension1",ylab="Scaling Dimension2",col = c(3,2)[g1],pch = 16)

######### Multiple comparison ###########
biocLite("qvalue"); library(qvalue);
q_obj <- qvalue(genediff[,"p_value"], fdr.level = 0.05) ## q value
hist(q_obj$pvalues,main = "", xlab = "p-value"); abline(v = 0.05/nrow(genediff),col = 2)
hist(q_obj$qvalues,main = "", xlab = "q-value"); abline(v = 0.05,col = 2)

length(which(q_obj$pvalues < 0.05/nrow(genediff))) ## FWER Differentially Expressed Gene
length(which(q_obj$qvalues < 0.05)) # FDR Differentially Expressed Gene

biocLite("hgu133plus2.db"); library( hgu133plus2.db ) ## gene annotation + explanation
genesymbol <- data.frame(symbol = sapply(contents(hgu133plus2SYMBOL), paste),
description = sapply(contents(hgu133plus2GENENAME),paste) )

rownames(genediff)[which(q_obj$qvalues < 0.05)] ## FDR DEG ID
genesymbol[which(q_obj$qvalues < 0.05),] ## Name of the gene

## END ## 