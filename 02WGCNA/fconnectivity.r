load("dataExpr.RData")
datExpr <- dataExpr
adjacency <- adjacency(datExpr, power = 1
moduleColors <- labels2colors(net$colors)
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
moduleConnectivity <- intramodularConnectivity(adjacency, moduleColors)
module <- "salmon"
moduleGenes <- which(moduleColors == module)
moduleConnectivitySubset <- moduleConnectivity[moduleGenes, ]

# 按连通性排序并选择hub基因
sortedConnectivity <- moduleConnectivitySubset[order(-moduleConnectivitySubset$kWithin), ]
write.table(sortedConnectivity,"module13_hubgenes.txt",sep="\t",quote=F,row.names=T,col.names=T)
