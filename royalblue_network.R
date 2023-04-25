library(dplyr)
load("./rdata/clustering/clust_sh.RData")
load("./rdata/datainput/annotation2016.RData")


dissTOM = as.data.frame(clust_sh$disTom)
adj = clust_sh$adj

colnames(adjacency) = clust_sh$probes
rownames(adjacency) = clust_sh$probes

colnames(dissTOM) = clust_sh$probes
rownames(dissTOM) = clust_sh$probes

royalblue = as.data.frame(adj[clust_sh$dynamicColors=='royalblue',clust_sh$dynamicColors=='royalblue'])
colnames(royalblue) = clust_sh$probes[clust_sh$dynamicColors=='royalblue']
rownames(royalblue) = clust_sh$probes[clust_sh$dynamicColors=='royalblue']

royalblue_probes = data.frame(ProbeName = (clust_sh$probes[clust_sh$dynamicColors=='royalblue']))
royalblue_probes_annot = inner_join(royalblue_probes, annotation)



diag(royalblue) <- NA
# quantile
q99 = quantile(royalblue, na.rm = TRUE, probs = 0.99)
q975 = quantile(royalblue, na.rm = TRUE, probs = 0.975)
q95 = quantile(royalblue, na.rm = TRUE, probs = 0.95)
q90 = quantile(royalblue, na.rm = TRUE, probs = 0.90)

# supprimer les valeurs inférieures à q99
royalblue_q99 = replace(royalblue, royalblue< q99, NA)
royalblue_q975 = replace(royalblue, royalblue< q975, NA)
royalblue_q95 = replace(royalblue, royalblue< q95, NA)
royalblue_q90 = replace(royalblue, royalblue< q90, NA)

#write.csv(royalblue_q99, "2023_03_31_royalblueq99.csv")
#write.csv(royalblue_q975, "2023_03_31_royalblueq975.csv")
#write.csv(royalblue_q95, "2023_03_31_royalblueq95.csv")
#write.csv(royalblue_q90, "2023_03_31_royalblueq90.csv")



royalblue_q95_bis = royalblue

for (i in 1:nrow(royalblue)) {
  s = (max(royalblue[i,], na.rm = TRUE))
 
  for (j in 1:ncol(royalblue)) {
    royalblue_q95_bis[i,j] =   replace(royalblue[i,j], royalblue[i,j]< q95, NA)  
    if (!is.na(royalblue[i, j]) && (s == royalblue[i, j]) && is.na(royalblue_q95_bis[i,j])) {
      royalblue_q95_bis[i,j] = royalblue[i, j] 
    }
  }
}

for (j in 1:ncol(royalblue)) {
  p = (max(royalblue[,j], na.rm = TRUE))
  for (i in 1:nrow(royalblue)) {
    if (!is.na(royalblue[i, j]) && (p == royalblue[i, j]) && is.na(royalblue_q95_bis[i,j])) {
      royalblue_q95_bis[i,j] = royalblue[i, j] 
      #print(p)
    }
  }
}


write.csv(royalblue_q95_bis, file = "2023_04_17_royalblue_q95.csv")

