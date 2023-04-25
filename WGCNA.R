######################################################
######################################################
###
###     PACKAGES INSTALLATION AND LOADING
###
######################################################
######################################################

# BiocManager::install("BiocInstaller")
# install.packages("BiocManager")
if (!require("WGCNA")) BiocManager::install("WGCNA")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()

######################################################
######################################################
###
###     CODE
###
######################################################
######################################################

######################################################
###     Soft Thresholding power
######################################################

PlotSoftThresholdingPower = function(data, networkType){
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))

  # Call the network topology analysis function
  sft = pickSoftThreshold(data, powerVector = powers, verbose = 5, networkType = networkType)

  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;

  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");

  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")

  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",
       ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

######################################################
###     Clustering
######################################################

Clustering.WGCNA = function(dataExpr = NULL, TOMRData = NULL, TOM = NULL, networkType = "unsigned", minModuleSize = NULL, softPower = 6, deepSplit = 2 ) {

  if (is.null(TOMRData) && is.null(TOM)){
    adjacency = adjacency(dataExpr, power = softPower, type = networkType);
    # Clustering using TOM
    # Turn adjacency into topological overlap
    TOM = TOMsimilarity(adjacency);
  } else if (!is.null(TOMRData)) {
    load(TOMRData)
  }
  dissTOM = 1-TOM
  probes = colnames(dataExpr)
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");

  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = deepSplit, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
  dynamicColors = labels2colors(dynamicMods);
  # Calculate eigengenes
  MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");

    return(list(probes = probes, geneTree=geneTree, dynamicMods = dynamicMods, dynamicColors = dynamicColors, MEs = MEs, METree = METree, dissTom =dissTOM, adj = adjacency))
}

######################################################
###
###     Merge Close Modules
###
######################################################

MergeCloseModules = function(data, dynamicColors, cutHeight = 0.25){
  MEDissThres = cutHeight
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(data, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
}

######################################################
###
###     Module-trait relationships
###
######################################################
ModuleTraitRelationships = function(datTraits, wgcna_data, cat = ""){
  nGenes = length(wgcna_data$probes);
  nSamples = nrow(datTraits);

  MEs = orderMEs(wgcna_data$MEs)

  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(5, 9, 3, 3));

  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships\n", cat, sep=' '))
}

######################################################
###
###     Module composition
###
######################################################

ModuleComposition.Probes = function(clustering, cat){
  modulesNames = levels(factor(clustering$dynamicColors))
  modulesList_probes =  vector(mode = "list", length = length(modulesNames))
  names(modulesList_probes) = paste(cat,modulesNames)

  idx = 0
  for (module in modulesNames){
    idx = idx+1
    temp1 = clustering$probes[which(clustering$dynamicColors==module)]
    modulesList_probes[idx] = list(temp1)

  }
  return (modulesList_probes);
}

ModuleComposition.Genes = function(clustering, cat, annotation){
  modulesNames = levels(factor(clustering$dynamicColors))
  modulesList =  vector(mode = "list", length = length(modulesNames))
  modulesList_genes = vector(mode = "list", length = length(modulesNames))
  names(modulesList_genes) = paste(cat,modulesNames)
  names(modulesList) = modulesNames

  idx = 0
  for (module in modulesNames){
    idx = idx+1
    temp = as.data.frame(clustering$probes[which(clustering$dynamicColors==module)])
    colnames(temp) = "ProbeName"
    temp = inner_join(temp,annotation)
    modulesList[idx] = list(temp,module)
    modulesList_genes[idx] = list(temp$Annotation)
  }
  return (modulesList_genes);
}
