# Match two columns of different dataframes and save in a new dataframe
setwd("E:/ALINE ARTIGO EPILEPSIA")


# To match the collumns - genes and miRNA
#mirtarbase
epilepsymirnamirtarbase <- merge(mirtarbase, genes,
             by.x = "Target Gene", by.y = "V1" )

#To save the dataframe in the cytoscape format
write.table(epilepsymirnamirtarbase,file="epilepsymirnamirtarbase.tsv",sep="\t",row.names = FALSE,quote=FALSE)

#tarbase
epilepsymirnatarbase <- merge(tarbase, genes,
                            by.x = "geneName", by.y = "V1" )

#To save the dataframe in the cytoscape format
write.table(epilepsymirnatarbase,file="epilepsymirnatarbase.tsv",sep="\t",row.names = FALSE,quote=FALSE)

