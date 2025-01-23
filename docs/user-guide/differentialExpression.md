# Running Differential Expression

_Adapted from the [Seurat differential expression vignette](https://satijalab.org/seurat/articles/de_vignette)_

When running differential expression on the Seurat object produced by SINCLAIR, these steps should be followed to ensure an accurate output

## Choosing the identities

Before running any sort of differential expression, the identities of the Seurat object need to be determined and set by the user. For example, the following code will set the Seurat object identities to the clusters found at resolution 0.8:

```
Idents(seuratObject) = "SCT_snn_res.0.8"
```

The identities are not limited to the clusters; these can be set to any categorical variable in the metadata, including cell types that have been determined by SingleR or classified by the user, cell cycle phase (G1/G2M/S), or experimental group.

## Preparing the Seurat object for differential expression with `PrepSCTFindMarkers`

After defining the object identities, Seurat requires that when when running differential expression on the SCT assay, the object needs to be prepared with the `PrepSCTFindMarkers` function:

```
seuratObject = PrepSCTFindMarkers(seuratObject)
```

The purpose of the `PrepSCTFindMarkers` is to use the minimum of the median UMI of individual objects to appropriately scale the SCT assay prior to differential expression<sup>[ref](#refPrepSCTfindMarkers)</sup>.

## Finding differentially expressed genes in Seurat

Two options are available when running differential expression through Seurat.

### Option 1: Running `FindMarkers`<sup>[ref](#refFindMarkers)</sup>

The "traditional" method of running differential expression compares two separate categories; this can still be exercised in the presence of more than the two comparison groups of interest. Following the example above, comparing the two largest clusters as determined by Seurat would be run as follows:

```
deGeneList = FindMarkers(seuratObject,ident.1=0,ident.2=1,test.use = "MAST")
```

By default, Seurat uses the non-parametric Wilcoxon rank-sum test to identify significant genes. The MAST<sup>[ref](#refMAST)</sup> algorithm uses a hurdle model to account for the sparsity of scRNASeq count matrices, and tends to be more sensitive to likely changes.

### Option 2: Running `FindAllMarkers`<sup>[ref](#refFindAllMarkers)</sup>

When looking to define a set of potential clusters or identities using gene markers, the `FindAllMarkers` function will run the differential expression by setting the first identity group (equivalent to `ident.1`) to each cluster in turn and using the remaining clusters as the comparison group (`ident.2`).

```
markerGeneList = FindAllMarkers(seuratObject,test.use="MAST")
```

The resulting list will show the differentially expressed genes that are significantly enriched or depleted for each identity within the preselected category. In essence, `FindAllMarkers` runs a loop where `FindMarkers` is run for each individual identity.

### Alternative methods for differential expression

The default values for some of the more frequently altered parameters for the `FindMarkers` and `FindAllMarkers` functions are as follows:

```
deGenes = FindMarkers(seuratObject, ident.1=NA, ident.2=NA, features=NULL,
  test.use="wilcox," logfc.threshold=0.25, min.pct=0.1)

```

If `ident.1` is defined and `ident.2` is left as a NULL value, the cells not included in `ident.1` will be used as the comparison, i.e. essentially running a one vs. rest comparison, similar to what is described in `FindAllMarkers`.

The `features` parameter can take a vector of specified genes and only run differential expression for those genes. Note that this assumes that the gene is not filtered out by any other criteria, such as those listed below. In all other cases, all genes are run through initial filtering and statistical testing.

The `test.use` parameter is used to change the statistical test applied to identify differentially expressed genes. As mentioned above, the MAST algorithm is often applied to single cell data, while other tests that can be applied include `"bimod"`, `"roc"`, `"t"`, `"negbinom"`, `"poisson"`, `"LR"`, and `"DESeq2"`.

The `logfc.threshold` filters out all genes that do not meet a log2 fold change threshold. This is applied to improve the speed of the computation, as those genes that do not meet this value are not subjected to statistical testing. The default threshold value is 0.25, which is roughly a 20% average fold change.

The `min.pct` threshold filters out genes that are not expressed in the fraction of cells below the threshold for **both** comparison groups. If only one comparison group is below the threshold, the statistics for the gene will still be calculated. The default `min.pct` threshold is 0.1, or 10%.

Since a number of genes will be unaffected by the experiment, a full differential expression list will require setting both `logfc.threshold` and `min.pct` to 0, so as to calculate the statistics for less relevant genes. These parameters may also need to be adjusted when selecting specific genes with the `features` parameter, as any user-defined genes that do not meet the threshold criteria will not be submitted to statistical testing.

## Expected outputs

The output of the `deGeneList` variable above will have a table structure resembling the following:

|       | p_val                | avg_log2FC        | pct.1 | pct.2 | p_val_adj            |
| ----- | -------------------- | ----------------- | ----- | ----- | -------------------- |
| Tff1  | 2.27649774378448e-20 | 0.503646816090369 | 0.313 | 0.039 | 5.15945448651315e-16 |
| Gkn1  | 2.37538488060338e-18 | 0.31537227402121  | 0.23  | 0.007 | 5.3835722933995e-14  |
| Gkn2  | 9.91735374467266e-14 | 0.238583210040284 | 0.166 | 0.004 | 2.24766905269261e-09 |
| Oaz1  | 1.26240861948424e-12 | 0.198987799763473 | 0.996 | 0.968 | 2.86112289519909e-08 |
| Rab5c | 3.39748073842636e-12 | 0.301705843310934 | 0.645 | 0.4   | 7.7000503455695e-08  |
| Rbm3  | 5.53653902770416e-12 | 0.226144933696604 | 0.97  | 0.893 | 1.25480120523887e-07 |
| Cfl1  | 9.22315807122862e-12 | 0.171800659275371 | 1     | 0.982 | 2.09033654526325e-07 |
| Cd83  | 2.27101797312635e-10 | 0.242223105538207 | 0.849 | 0.789 | 5.14703513429356e-06 |
| Tff2  | 2.66531282660886e-10 | 0.226607061165825 | 0.189 | 0.025 | 6.04066499022633e-06 |

The first unlabeled column lists the genes, followed by raw p-value, average log2 fold change, the percentage of cells expressing the gene in populations 1 and 2, and finally the adjusted p-value. The sign of the fold change and the definitions of `pct.1` and `pct.2` are defined by the order of groups selected in the `FindMarkers` call: Positive fold changes indicate enrichment in the first group (i.e. defined as `ident.1`), as well as the `pct.1` value

## Common issues and questions

#### Why are there so many genes with extremely small p-values?

In most circumstances, a small p-value generally indicates that the gene is extremely significant and should not be ignored. However, many of the statistical tests that are designed based on distributions, including MAST, operate under a set of assumptions. With single cell, one of the assumptions that gets violated is an upper limit on the number of replicates included, since each individual cell is treated as a replicate and leads to a comparison of 2 groups with as many as tens of thousands of replicates each.

Unfortunately, it falls to the user to determine the relative significance of the genes that are identified as significant by examining the other statistics provided (i.e. `avg_log2FC` and `pct.1`/`pct.2`). Additionally, genes of interest can also be run through the `AverageExpression` and the `FeaturePlot` functions to explore the overall expression of the gene in the individual contrast groups. Other approaches, such as pseudobulk differential expression (see below) might be implemented in order to address this p-value phenomenon

#### Requiring `JoinLayers`

Since version 5, Seurat keeps individual samples as "layers" in the S4 data structure. This makes it simpler to apply various functions to each sample, such as SCTransform normalization or variable feature identification, since a single call to the Seurat object will behave like the `lapply` function in R. However, this also keeps the individual counts tables separate, which makes differential expression nigh impossible. To address this, the user needs to join the layers prior to running differential expression:

```
so_dePrep = JoinLayers(seuratObject)
so_dePrep = PrepSCTFindMarkers(so_dePrep)
markers = FindAllMarkers(so_dePrep,test.use="MAST")
```

#### Running differential expression on subsets of cells

Oftentimes the user will be interested in examining two different subsets of cells; for example, comparing the two experimental contrast groups in B cells alone. When isolating a group with the subset in question, the user might send a command such as:

```
Bcells = subset(seuratObject,cells=colnames(seuratObject)[which(so$cellType=="B cells")])
```

However, when running `FindMarkers` on this subset, the following may be encountered:

```
FindMarkers(Bcells,ident.1="group1",test.use="MAST")

Error in FindMarkers.SCTAssay(object = data.use, slot = slot, cells.1 = cells$cells.1,  :
  Object contains multiple models with unequal library sizes. Run `PrepSCTFindMarkers()` before running `FindMarkers()`.

```

In this event, the user should ensure that the subset was extracted from the Seurat object that has already been prepared through `PrepSCTFindMarkers` to ensure appropriate scaling. The second step is to use a hidden parameter `recorrect_umi` to tell the program to "ignore" the scaling step.

```
FindMarkers(Bcells,ident.1="group1",test.use="MAST", recorrect_umi=FALSE)
```

In theory, this flag can also be used in all other steps in order to skip the scaling step for differential expression, but this is generally not recommended.

#### Pseudobulk differential expression

If there are enough individual replicates, a user can convert the single cell dataset into a pseudobulk dataset and treat the individual samples as pooled replicates to be submitted through a standard RNASeq differential expression protocol, such as Limma or DESeq2. This has been explored using the [Libra](https://github.com/neurorestore/Libra) tool and through [Seurat](https://satijalab.org/seurat/articles/de_vignette#perform-de-analysis-after-pseudobulking) using their `AggregateExpression` function. The results tend to increase the p-values, which subsequently reduces the likelihood of false positives in identifying differentially expressed genes.

The main warning when running a pseudobulk approach is to ensure that there are enough samples to warrant aggregation; if only one sample is available per experimental condition, the user will be attempting a 1v1 differential expression design, which is not nearly robust enough to account for any intra-group sample variability.

</br>
</br>

<font size='2'>References:

<a name=refPrepSCTfindMarkers>1.</a> [PrepSCTFindMarkers](https://satijalab.org/seurat/reference/prepsctfindmarkers)

<a name=refFindMarkers>2.</a> [FindMarkers](https://satijalab.org/seurat/reference/findmarkers)

<a name=refFindAllMarkers>3.</a> [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers)

</font>

<font size='2'>Author: Nathan Wong. January 2023</font>
