## Tutorial retrieved from https://www.mirnet.ca/miRNet/docs/RTutorial.xhtml jan2022

# Step 0: install dependencies

install.packages("pacman")

library(pacman)

pacman::p_load(RSQLite, Cairo, fastmatch, igraph, RJSONIO, foreach, doParallel, preprocessCore, limma, edgeR, HTqPCR, genefilter)

# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# Step 2: Install miRNetR WITHOUT documentation
devtools::install_github("xia-lab/miRNetR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))

# Step 2: Install miRNetR WITH documentation
devtools::install_github("xia-lab/miRNetR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

library(miRNetR)

#### Step 1. Initiate the dataSet object
Init.Data("mir", "mirlist")
#> "miRNetR init done!"

#### Step 2. Set up the user input data
SetupMirListData(mirs = "hsa-mir-101-3p
hsa-mir-133b
hsa-mir-147a
hsa-mir-3140-3p
hsa-mir-361-5p
hsa-mir-510-5p", orgType = "hsa", idType = "mir_id", tissue = "Kidney")

#### Step 3. Set up targets
nms.vec = c("gene")
SetCurrentDataMulti()

#### Step 4. Perform miRNAs to target genes mapping, 
#### results are downloaded in your working directory ("mirnet_mir_target.csv" and "mir_tissue.csv")
QueryMultiListMir()
head(dataSet$mir.res, n = 3L)


#### Step 5. Generate miRNA-gene network files
CreateMirNets(net.type = "mir2gene")

#### Step 6. Prepare network files, 
#### results are downloaded in your working directory ("node_table_mirnet_0.csv", "mirnet_0.json" and "mirnet.graphml")
PrepareMirNet(mir.nm = "mirnet1", file.nm = "mirnet_0.json")

#### Step 7. Perform miRNA family enrichment analysis, 
#### results are downloaded in your working directory ("network_enrichment_mirfamily_1.json" and "mirnet_enrichment.csv")
PerformMirTargetEnrichAnalysis(
  adjust.type = "NA",
  fun.type = "mirfamily",
  file.nm = "network_enrichment_mirfamily_1",
  IDs = "hsa-mir-101-3p; hsa-mir-147a; hsa-mir-361-5p; hsa-mir-133b; hsa-mir-510-5p; hsa-mir-3140-3p",
  algo = "hyp"
)

resTable <- read.csv("mirnet_enrichment.csv", header=T, as.is=T)
head(resTable, n = 3L)

###### Mapping to multiple targets



