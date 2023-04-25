# Small networks of expressed genes in the whole blood and relationships to profiles in circulating metabolites provide insights in inter-individual variability of feed efficiency in growing pigs
Camille Juigné, Emmanuelle Becker and Florence Gondret

*RMarkdown notebooks*

1) Import and process data : [input_data_processing](https://github.com/cjuigne/multiomics_and_feed_efficiency/blob/main/input_data_processing.Rmd)
1) (1-bis) RData files for input data (loaded in the next notebooks) [./rdata/datainput](https://github.com/cjuigne/multiomics_and_feed_efficiency/tree/main/rdata/datainput)
2) Then, we constructed the weighted gene co-expression networks using WGCNA : [wgcna](https://github.com/cjuigne/multiomics_and_feed_efficiency/blob/main/wgcna.Rmd)
3) We established profiles of circulating metabolites: [metabolites](https://github.com/cjuigne/multiomics_and_feed_efficiency/blob/main/pca_met.Rmd), [fatty acids](https://github.com/cjuigne/multiomics_and_feed_efficiency/blob/main/pca_correlation_wgcna.Rmd)
4) and evaluating connections between metabolic and transcriptomic levels: [pca_correlation_wgcna](https://github.com/cjuigne/multiomics_and_feed_efficiency/blob/main/pca_correlation_wgcna.Rmd)
