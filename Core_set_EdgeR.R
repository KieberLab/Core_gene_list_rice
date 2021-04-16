#Rice BA treated, RNAseq.

setwd()
getwd()
getwd()
library(edgeR)
library(limma)
library(reshape2)
library(grid)
library(stringr)
library(ggplot2)
library(GenomeInfoDb)
library(DESeq2) #need to install GenomeInfoDbData
library(GGally)


x=read.delim("new_Potter_2_counts.txt")
head(x)
tail(x)


#reading counts for a specific gene

x[x[,1]=="LOC_Os01g01010",]

#adding a grouping factor

group=factor(c(1,1,1,2,2,2,3,3,3))

?DGEList
y <- DGEList(counts=x[,2:10], genes=x[,1], group=group)

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
y=calcNormFactors(y, method="TMM")
keep=rowSums(cpm(y)>2)>=3 #keeping genes w/ reasonable expression (higher than 2cpm in at least 3 libraries )
y=y[keep, ]


#checking if the filtering worked
dim(x)
dim(y)

#this checks counts of a specific gene

x[x[,1]=="LOC_Os01g09260",]
x[x[,1]=="LOC_Os01g01090",]

y.cor.cpm[y.cor.cpm[,1]=="LOC_Os01g09260",]

#this is an important step in generating the glm. I used the codes provided in the edgeR manual
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
y.cor.cpm=data.frame(y$genes, cpm(y))
head(y.cor.cpm)
y.cor.cpm[y.cor.cpm[,1]=="LOC_Os01g09260",]

#summarize it in a table
write.table(y.cor.cpm, "y.cor.cpm_potter2.txt", sep="\t", row.names = F, quote = F)
#looking at the CPM of specific gene
y.cor.cpm[y.cor.cpm[,1]=="LOC_Os01g01080",]

#this is an MDS plot. The colors 1:8 will be repeated three times (3 for the control and 3 for BA (I want them to be the same color). We have 48 samples so 48/6 is 8, therefore 1:8 in the argument "col". For the point(plotting character "pch") - I want two different point, for example circles (16) and triangles (17). Therefore, I will have 16:17, each 3 times, namely 16 16 16 (circle, circle, circle) and then 17 17 17 (triangle triangle triangle). R repeats it until it's done with the whole string. The legend is a little tricky - in this case each color will be repeated 2 times (I don't want to have all the data point on my legend!). Same logic goes for the points 16:17 each 1 times - 16 17 16 17 16 17...and so on (circle triangle circle triangle....)
?plotMDS

plotMDS(y, col=rep(1:3, each=3), labels=NULL, pch=rep(16:18, each=3), legend=TRUE , legend(("bottomleft"),legend=c("0h_control", "2h_control", "2h_BA") , col=rep(1:3, each=1), pch=rep(16:18, each=1, cex=3)))



#DE analysis.Now it's time to generate the tables with the time point comparisons

#glm for all the time-points
?glmLRT
lrt.2hBA_2hC <- glmLRT(fit, contrast=c(0,-1,1)) #BA_2h vs control_2h
lrt.120C_0hC <- glmLRT(fit, contrast=c(-1,1,0)) #control_2h vs control_0h
lrt.120BA_0hC <- glmLRT(fit, contrast=c(-1,0,1)) #BA_2h vs control_0h

#"toptags"
?topTags

tt.2hBA_2hC <-topTags(lrt.2hBA_2hC, n=20278)
tt.120C_0hC <-topTags(lrt.120C_0hC, n=20278)
tt.120BA_0hC <-topTags(lrt.120BA_0hC, n=20278)



write.table(tt.2hBA_2hC$table, "tt.2hBA_2hC_potter2.txt", sep="\t", row.names = F, quote = F)
write.table(tt.120C_0hC$table, "tt.2hC_0hC_potter2.txt", sep="\t", row.names = F, quote = F)
write.table(tt.120BA_0hC$table, "tt.2hBA_0hC_potter2.txt", sep="\t", row.names = F, quote = F)


#selecting w/ logFC threshold
logFC=0.7
FDR=0.05

###positive logFC###
BA2h_2hC_2h_UP=tt.2hBA_2hC[which(tt.2hBA_2hC$table$logFC > logFC & tt.2hBA_2hC$table$FDR < FDR),] #... genes
write.table(BA2h_2hC_2h_UP$table, "1_Potter2_2h_UP.txt", sep="\t", row.names=F, quote=F)
dim(BA2h_2hC_2h_UP)
C120_0hC_UP=tt.120C_0hC[which(tt.120C_0hC$table$logFC > logFC & tt.120C_0hC$table$FDR < FDR),] #... genes
write.table(C120_0hC_UP$table, "Potter2_time_control_UP.txt", sep="\t", row.names=F, quote=F)
BA120_0hC_UP=tt.120BA_0hC[which(tt.120BA_0hC$table$logFC > logFC & tt.120BA_0hC$table$FDR < FDR),] #... genes
write.table(BA120_0hC_UP$table, "Potter2_time_BA_UP.txt", sep="\t", row.names=F, quote=F)

###negartive logFC###
BA2h_2hC_2h_down=tt.2hBA_2hC[which(tt.2hBA_2hC$table$logFC < -logFC & tt.2hBA_2hC$table$FDR < FDR),] #... genes
write.table(BA2h_2hC_2h_down$table, "1_Potter2_2h_down.txt", sep="\t", row.names=F, quote=F)
dim(BA2h_2hC_2h_down)

C120_0hC_down=tt.120C_0hC[which(tt.120C_0hC$table$logFC < -logFC & tt.120C_0hC$table$FDR < FDR),] 
write.table(C120_0hC_down$table, "Potter2_time_control_down.txt", sep="\t", row.names=F, quote=F)
BA120_0hC_down=tt.120BA_0hC[which(tt.120BA_0hC$table$logFC < -logFC & tt.120BA_0hC$table$FDR < FDR),] 
write.table(BA120_0hC_down$table, "Potter2_time_BA_down.txt", sep="\t", row.names=F, quote=F)

dim(BA2h_2hC_2h_UP)
dim(C120_0hC_UP)
dim(BA120_0hC_UP)
dim(BA2h_2hC_2h_down)
dim(C120_0hC_down)
dim(BA120_0hC_down)



