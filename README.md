# nsmblR: Network inference ensemble modeling

nsmblR is an algorithm that infers a consensus gene regulatory network based on a gene expression data compendium. This tool uses a total of 7 different inference algorithms and performs a voting of the edges discovered by each algorithm to achieve a consensus voting on the relevant edges. The algorithms currently used by nsmblR are:

  - [CLR](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050008): Mutual information based algorithm with distribution correction
  - [ARACNE](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7): Mutual information based algorithm with edge pruning
  - Spearman correlation
  - [PCIT](https://www.ncbi.nlm.nih.gov/pubmed/20007253): Partial correlation based algorithm
 - [MRNET](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3171353/): Maximum relevance/minimum redundancy network inference based on mutual information
 - [MRNETB](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3171353/)
 - [MutRank](https://www.ncbi.nlm.nih.gov/pubmed/19767600): Rank correlation inference
 
Please follow the links above for the appropriate references and algorithm descriptions. After running all of the algorithms, the inferred results are filtered based on the quantile of the scores (default threshold set to 0.97) and the edges are tallied. Multiple voting schemes are provided and the final results account for edges that are present in more than 51% of the cases (4 algorithms out of 7.) 

## Installation

Install from GitHub using `devtools` as:

```
devtools::install_github("diogocamacho/nsmblR")
```

## Running
The easiest way to run `nsmblR` is to use its wrapper `ensemble_model` as:

```
library(nsmblR)
res <- ensmeble_model(data, gene_names)
```

where `data` is a gene expression compendium (genes on rows, samples on columns) and `gene_names` are the gene symbols for the rows. Internally `nsmblR` will subset the final set of edges as those that are consistent across multiple inference methods that are also in the top 97% quantile of the edge score. For greater flexibility, look into the individual functions of the package to taylor the results to your specific problem.