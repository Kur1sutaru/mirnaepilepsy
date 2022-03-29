library("edgeR")

setwd("E:/ALINE ARTIGO EPILEPSIA/GSE127871")


files <- c("Sample_31.txt",
           "Sample_38.txt",
           "Sample_55.txt",
           "Sample_56.txt",
           "Sample_60.txt",
           "Sample_63.txt",
           "Sample_65.txt",
           "Sample_67.txt",
           "Sample_68.txt",
           "Sample_80.txt",
           "Sample_81.txt",
           "Sample_82.txt")

DG <- readDGE(files,header=TRUE)

#To create a data frame now

#To load each sample individually

sample_31 = read.table('Sample_31.txt', sep = '\t',
                   header = T)
sample_38 = read.table('Sample_38.txt', sep = '\t',
                   header = T)
sample_55 = read.table('Sample_55.txt', sep = '\t',
                           header = T)
sample_56 = read.table('Sample_56.txt', sep = '\t',
                            header = T)
sample_60 = read.table('Sample_60.txt', sep = '\t',
                            header = T)
sample_63 = read.table('Sample_63.txt', sep = '\t',
                            header = T)
sample_65 = read.table('Sample_65.txt', sep = '\t',
                          header = T)
sample_67 = read.table('Sample_67.txt', sep = '\t',
                            header = T)
sample_68 = read.table('Sample_68.txt', sep = '\t',
                            header = T)
sample_80 = read.table('Sample_80.txt', sep = '\t',
                            header = T)
sample_81 = read.table('Sample_81.txt', sep = '\t',
                       header = T)
sample_82 = read.table('Sample_82.txt', sep = '\t',
                       header = T)

# To combine data sets into a matrix
geneCounts = data.frame(sample_31 [,2], sample_38 [,2],
                        sample_55 [,2], sample_56 [,2],
                        sample_60 [,2], sample_63 [,2],
                        sample_65 [,2], sample_67 [,2],
                        sample_68 [,2], sample_80 [,2],
                        sample_81 [,2], sample_82 [,2])
           
row.names(geneCounts) = sample_80[,1]
sizeGeneCounts = dim(geneCounts)
geneCounts = geneCounts[1:(sizeGeneCounts[1]-5),]
condition = c(rep('low',6), rep('high',6))
sampleNames = c('sample_31', 'sample_38', 'sample_55',
                'sample_56', 'sample_60','sample_63', 
                'sample_65', 'sample_67', 'sample_68',
                'sample_80', 'sample_81', 'sample_82')
colnames(geneCounts) = sampleNames
View(geneCounts)

#To build the generalized linear model that will be used for differential expression 
dge <- DGEList(counts=geneCounts, group=condition)
design <- model.matrix(~condition+0, data=dge$samples)
colnames(design) = gsub("condition","",colnames(design))

# To perform the normalization 
dge <- calcNormFactors(dge)
dge$samples
plotMDS(dge)

#To write a plot (PCA)
jpeg(file="GSE127871.jpeg", width=5000, height=5000, units="px", res=300)
plotMDS(dge)
dev.off()

#To estimate the dispersion
disp <- estimateGLMCommonDisp(dge, design)
disp <- estimateGLMTrendedDisp(disp, design)
disp <- estimateGLMTagwiseDisp(disp, design)
plotBCV(disp)

#To analyze the differential expression and likelihood ratio
fit <- glmQLFit(disp,design)
qlf <- glmQLFTest(fit,coef=2)

fit <- glmFit(disp, design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

#To confirm if conditions are correct
colnames(design)

#To perform the exact test
dgeTest <- exactTest(disp)
dgeTest

#To see the histogram of P-Values
hist(dgeTest$table[,"PValue"], breaks=50)


#To tell edgeR what comparison you want to perform
lowxhigh = makeContrasts(low-high, levels=design)


##### seizures/month low x seizures/month high ####


#To perform the differential expression testing for that comparison
lrt.lowxhigh = glmLRT(fit, contrast=lowxhigh)
res_lowxhigh<-topTags(lrt.lowxhigh, n=500000, sort.by = "p.value")
write.csv(res_lowxhigh, file="lowxhigh_GSE127871.csv")

#To write a table with the significative p-values
sig_lowxhigh<-topTags(lrt.lowxhigh, n=500000, p.value=.05, sort.by = "p.value")
write.csv(sig_lowxhigh, file="lowxhigh.csv")
