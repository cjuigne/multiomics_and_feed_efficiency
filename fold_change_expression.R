if (!require("ggplot2")) install.packages("ggplot2")

FoldChangeExpr = function(dataExpr, threshold){
  data = data.frame(rownames(dataExpr), row.names = rownames(dataExpr))
  data$min = apply(t(dataExpr), 2, min)
  data$max = apply(t(dataExpr), 2, max)
  data$FC =  data$max - data$min 
  select = which(data$FC>threshold)
  return(dataExpr[select,])
}

FoldChangeExpr.distribution = function(dataExpr,cat =""){
  data = data.frame(rownames(dataExpr), row.names = rownames(dataExpr))
  data$min = apply(t(dataExpr), 2, min)
  data$max = apply(t(dataExpr), 2, max)
  data$FC =  data$max - data$min 
  return(ggplot(data, aes(FC)) + 
    geom_histogram(colour="black", fill="white", binwidth = 0.25) +geom_vline(xintercept=1, color = 'red') +
    labs(title=paste(" Distribution of gene expression",cat, sep = " ")))
}
