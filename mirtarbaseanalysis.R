setwd("C:/miRNAsDoc/Analise_miRNAsValidados/")

library(readxl)

##### Auxiliar functions ######

## filter mir-tg list using seed for Regulator
filterInteractionsUsingSeedForRegulator <- function(interactions,seed.reg) {
  ## all interaction data follows this default organization: 
  ## first column is regulator, second column is target
  
  keepRows <- which(is.element(interactions[,1],seed.reg))
  interactions <- interactions[keepRows,]
  interactions
}

standardizeGeneSymbols <- function(geneList,mapFunction,uppercase=FALSE){
  if(uppercase) {
    geneList <- toupper(geneList)
  }
  temp<-mapFunction[as.character(geneList)]
  temp <- as.character(temp)
  return(temp)
}

####### ----------------########



### Load datasets

### miRTarBase release 8
mirtarbase<- read_excel("mirtarbase_release8_hsa_MTI.xlsx") ##setup file path

#adjust column orders
mirtarbase<- as.data.frame(mirtarbase[,c("miRNA","Target Gene", "Target Gene (Entrez Gene ID)","miRTarBase ID","Support Type")])


### TarBase v8
tarbase<- read.delim("TarBase_v8_download.txt") ##setup file path
tarbase <- tarbase[which(tarbase$species == "Homo sapiens"),]
tarbase <- tarbase[,c("mirna","geneName","positive_negative","direct_indirect","up_down","method")]
tarbase$geneName <- gsub("\\(hsa\\)","",tarbase$geneName)

## HGCN mapping function (defined in Postdoc/Proheto_CardiacHypertrophy_SysBio)
load("HGCN_MapFunction.RData")



### Define miRNAs of interest
mir.interest <- c("hsa-miR-122-5p","hsa-miR-126-3p", "hsa-miR-150-5p", "hsa-miR-192-5p","hsa-miR-146a-5p")

##Query miRNAs targets from miRTarBase

targets.mirtarbase <- filterInteractionsUsingSeedForRegulator(mirtarbase,mir.interest)
targets.mirtarbase <- targets.mirtarbase[which(targets.mirtarbase$`Support Type`!="Non-Functional MTI"),]
targets.mirtarbase <- unique(targets.mirtarbase[,c("miRNA","Target Gene","Support Type")])
colnames(targets.mirtarbase) <- c("miRNA","gene","SupportType")                                         

##Query miRNAs targets from TarBase
targets.tarbase <- filterInteractionsUsingSeedForRegulator(tarbase,mir.interest)
targets.tarbase <- targets.tarbase[which(targets.tarbase$positive_negative=="POSITIVE"),]
targets.tarbase <- targets.tarbase[which(targets.tarbase$direct_indirect=="DIRECT"),]
targets.tarbase <- targets.tarbase[which(targets.tarbase$up_down=="DOWN"),]

## HITS-CLIP and PAR-CLIP are considered "limited evidence": keep only interactions supported by 2/4 or more experiments
## OBS.: Usually I keep those supported by >= evidences, nonetheless, in this case we have a miRNA  apparently more investigated than the
## other. A high threshold would result in few targets for one miRNA; a low threshold would result in an excessive 
## number of targets for the other. (See the comment below)
targets.tarbase.filt <- targets.tarbase[-which(targets.tarbase$method %in% c("HITS-CLIP","PAR-CLIP")),]
targets.tarbase.filt <- unique(targets.tarbase.filt[,c("mirna","geneName")])
targets.tarbase.weak <- targets.tarbase[which(targets.tarbase$method %in% c("HITS-CLIP","PAR-CLIP")),]
targets.tarbase.weak <- plyr::count(targets.tarbase.weak,c("mirna","geneName"))

## one miRNA has a much higher number of targets than the other (2676 vs 67), which may indicate some publication bias.
## thus, for miR-451a we keep all targets with >= 2 supporting evidence, whereas for miR-22-3p with >=4 supporting evidences
targets.tarbase.weak <- targets.tarbase.weak[which(targets.tarbase.weak$freq>=2),]
miR1225p <- plyr::count(targets.tarbase.weak[which(targets.tarbase.weak$mirna=="hsa-miR-122-5p"),],c("mirna","geneName"))
miR1263p <- plyr::count(targets.tarbase.weak[which(targets.tarbase.weak$mirna=="hsa-miR-126-3p"),],c("mirna","geneName"))
miR1505p <- plyr::count(targets.tarbase.weak[which(targets.tarbase.weak$mirna=="hsa-miR-150-5p"),],c("mirna","geneName"))
miR1925p <- plyr::count(targets.tarbase.weak[which(targets.tarbase.weak$mirna=="hsa-miR-192-5p"),],c("mirna","geneName"))
miR146a5p <- plyr::count(targets.tarbase.weak[which(targets.tarbase.weak$mirna=="hsa-miR-146a-5p"),],c("mirna","geneName"))

targets.tarbase.weak <- rbind(miR1225p[which(miR1225p$freq>=5),1:2],
                              miR1263p[which(miR1263p$freq>=5),1:2],
                              miR1505p[which(miR1505p$freq>=5),1:2],
                              miR1925p[which(miR1925p$freq>=5),1:2],
                              miR146a5p[which(miR146a5p$freq>=5),1:2])
## using the same annotation adopted by miRTarBase
targets.tarbase2 <- unique(rbind(cbind(targets.tarbase.filt[,1:2],SupportType="Functional MTI"),
                                 cbind(targets.tarbase.weak[,1:2],SupportType="Functional MTI (Weak)")))
colnames(targets.tarbase2) <- c("miRNA","gene","SupportType")
targets.tarbase <- targets.tarbase2
rm(targets.tarbase2)
rm(miR1225p,miR1263p,miR1505p,miR1925p,miR146a5p)

##merge datasets, removing redundancy
targets.mirs <- unique(rbind(targets.mirtarbase,targets.tarbase))

## standardize gene symbols using HGNC nomenclature
targets.mirs$gene <- standardizeGeneSymbols(targets.mirs$gene,genename.map,uppercase = FALSE)

## order targets by support type - strong support comes first
targets.mirs <- targets.mirs[order(targets.mirs$SupportType),]

##remove interactions supported by "Functional MTI" and "Functional MTI (Weak)", keeping only the first evidence
targets.mirs <- targets.mirs[-which(duplicated(targets.mirs[,1:2])),]

## number of targets per miRNA
table(targets.mirs$miRNA)


### export network and type of nodes - for Cytoscape
write.table(targets.mirs,file="mir_target_analysis_RS-OB.tsv",sep="\t",row.names = FALSE,quote=FALSE)

type.node <- rbind(cbind(unique(targets.mirs$miRNA),"miRNA"),
                   cbind(unique(targets.mirs$gene),"gene"))
write.table(type.node,file="mir_target_analysis_nodeType_RS-OB.tsv",sep="\t",row.names = FALSE,quote=FALSE)



### Run Functional Enrichment Analysis

## if needed, install clusterProfiler package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#
#BiocManager::install("clusterProfiler")

library("clusterProfiler")

## gene ids should be given as entrez id for functional enrichment functions
targets.mirs.entrez <- targets.mirs
eg <- bitr(targets.mirs$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id.map <- eg$ENTREZID
names(id.map) <- eg$SYMBOL
rm(eg)
targets.mirs.entrez$gene<- id.map[targets.mirs.entrez$gene]
targets.mirs.entrez<- targets.mirs.entrez[which(!is.na(targets.mirs.entrez$gene)),]


## create a list of targets for each miRNA - useful for functions that compare enrichment
## among groups of genes
list.targets <-sapply(mir.interest,function(x) unique(as.character(targets.mirs.entrez$gene[grep(x,targets.mirs.entrez$miRNA)])))

## run enrichment using KEGG pathway. compareCluster function provides interface to visually compare the enrichments
ck <- compareCluster(geneCluster = list.targets, fun = enrichKEGG)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

head(ck) 

dotplot(ck,showCategory=15, font.size=8)


write.table(ck@compareClusterResult,file="KEGG_Enrichment_Results_RS-OB.xls",sep="\t",row.names = FALSE)

save.image("findTargets_RS-OB.RData")