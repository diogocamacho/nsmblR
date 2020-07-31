# nsmblR: Network inference ensemble modeling

nsmblR is an algorithm that infers a consensus gene regulatory network based on a gene expression data compendium. This tool uses a total of 7 different inference algorithms and performs a voting of the edges discovered by each algorithm to achieve a consensus voting on the relevant edges. The algorithms currently used by nsmblR are:

  - [CLR](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050008): Mutual information based algorithm with distribution correction
  - [ARACNE](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7): Mutual information based algorithm with edge pruning
  - Spearman correlation
  - [MRNET](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3171353/): Maximum relevance/minimum redundancy network inference based on mutual information
  - [MRNETB](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3171353/)

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

where `data` is a gene expression compendium (genes on rows, samples on columns) and `gene_names` are the gene symbols for the rows. Internally `nsmblR` will subset the final set of edges as those that are consistent across multiple inference methods that are also in the top 97% quantile of the edge score. For greater flexibility, look into the individual functions of the package to tailor the results to your specific problem.

## Example

The following is an example on how to run the `nsmblR` package with data that is provided with the package. This data comes from _E. coli_ and it is a set of 200 genes in 20 different conditions. First, we will load the package:

```
library(nsmblR)
```

Now we will run the package wrapper on the example data.

```
net <- ensemble_model(data = data_matrix, gene_names = genes)
```

The `net` variable is a list that contains all inferred networks, the inferred networks filtered based on a quantile threshold (see documentation for the `edge_filtering` function), and the results on the consensus voting. The consensus vote is ultimately the result of the `nsmblR` package and can be access as:

```
consensus_net <- net$consensus_network
```

which is a data frame with N edges and 6 columns, where the columns are the gene names for the edges (x and y), and the presence of the edge in different voting regimens (majority, super_majority, absolute_majority, and quorum -- see documentation on the `edge_voting` for a detailed explanation of voting schemes.) 

You can then use a package like `igraph` to display the edges inferred. In this example, we can look at the edges that pass the super-majority condition as such:

```
library(igraph)
G <- graph_from_data_frame(N$consensus_network %>% dplyr::filter(., super_majority == 1), directed = FALSE)
plot(G)
```

which will show the edges that pass the super majority voting.

**NOTE**: this is only an example. Even as such, with a set of genes and samples selected at random, we see the inference of a relationship between the lsrR regulator and the lsrG gene, which is to be expected, thereby giving some confidence about the approaches used for inference and the voting methods employed here. (see more at [EcoCyc](https://ecocyc.org/gene?orgid=ECOLI&id=G6805))
