# GSEA API for Pandas
[![Build Status](https://travis-ci.com/krassowski/gsea-api.svg?branch=master)](https://travis-ci.com/krassowski/gsea-api) 
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](http://choosealicense.com/licenses/mit/)
[![DOI](https://zenodo.org/badge/188071398.svg)](https://zenodo.org/badge/latestdoi/188071398)

Pandas API for Gene Set Enrichment Analysis in Python (GSEApy, cudaGSEA, GSEA)

- aims to provide a unified API for various GSEA implementations; uses pandas DataFrames and a hierarchy of Pythonic classes.
- file exports (exporting input for GSEA) use low-level numpy functions and are much faster than in pandas
- aims to allow researchers to easily compare different implementations of GSEA, and to integrate those in projects which require high-performance GSEA (e.g. massive screening for drug-repositioning)
- provides useful utilities for work with GMT files, or gene sets and pathways in general in Python


## Installation

To install the API use:

```bash
pip3 install gsea_api
```

See [below](#Installing-GSEA-implementations) for the instructions on installation of specific GSEA implementations.

## Example usage

```python
from pandas import read_table
from gsea_api.expression_set import ExpressionSet
from gsea_api.gsea import GSEADesktop
from gsea_api.molecular_signatures_db import GeneSets

reactome_pathways = GeneSets.from_gmt('ReactomePathways.gmt')

gsea = GSEADesktop()

design = ['Disease', 'Disease', 'Disease', 'Control', 'Control', 'Control']
matrix = read_table('expression_data.tsv', index_col='Gene')

result = gsea.run(
    # note: contrast() is not necessary in this simple case
    ExpressionSet(matrix, design).contrast('Disease', 'Control'),
    reactome_pathways,
    metric='Signal2Noise',
    permutations=1000
)
```


Where `expression_data.tsv` is in the following format:

```
Gene	Patient_1	Patient_2	Patient_3	Patient_4	Patient_5	Patient_6
TACC2	0.2	0.1	0.4	0.6	0.7	2.1
TP53	2.3	0.2	2.1	2.0	0.3	0.6
```

### MSigDB integration

[Molecular Signatures Database](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) (MSigDB) can be downloaded from the [Broad Institute GSEA website](https://www.gsea-msigdb.org/gsea/downloads.jsp). It provides expert-curated gene set collections, as well as curated subset of pathway databases (Reactome, KEGG, Biocarta, Gene Ontology) trimmed to remove redundant, overlapping and and otherwise little-value terms (if needed).

You can download all the pathways collections at once (search for `ZIPped MSigDB` on the download page). After downloading and un-zipping (e.g. to a local directory named `msigdb`), you can access the gene sets from MSigDB with:

```python
from gsea_api.molecular_signatures_db import MolecularSignaturesDatabase

msigdb = MolecularSignaturesDatabase('msigdb', version=7.1)
msigdb.gene_sets
```

`msigdb.gene_sets` returns a list of dictionaries describing auto-detected pathways:

```python
[
    {'name': 'c1.all', 'id_type': 'symbols'},
    {'name': 'c1.all', 'id_type': 'entrez'},
    {'name': 'c2.cp.reactome', 'id_type': 'symbols'},
    {'name': 'c2.cp.reactome', 'id_type': 'entrez'}
    # etc..
]
```

Information about the location on disk and version are avilable in `msigdb.path` and `msigdb.version`.

`msigdb.load` loads the specific collection into a `GeneSets` object:

```python
> kegg_pathways = msigdb.load('c2.cp.kegg', 'symbols')
> print(kegg_pathways)
<GeneSets 'c2.cp.kegg' with 186 gene sets>
```

This object can be passed to any of the supporteed GSEA implementations; please see below for a detailed description of the `GeneSets` object.

### `GeneSets` objects

`GeneSets` represents a collection of sets of genes, where each set is represented as `GeneSet` object.

You can check the number of sets contained within a collection with:

```python
> len(kegg_pathways)
186
```

The gene sets are accessible with `gene_sets` (tuple) and `gene_sets_by_name` (dict) properties:

```python
> kegg_pathways.gene_sets[:2]
(<GeneSet 'KEGG_TIGHT_JUNCTION' with 132 genes>, <GeneSet 'KEGG_RNA_DEGRADATION' with 59 genes>)
> kegg_pathways.gene_sets_by_name
{
    'KEGG_TIGHT_JUNCTION': <GeneSet 'KEGG_TIGHT_JUNCTION' with 132 genes>,
    'KEGG_RNA_DEGRADATION': <GeneSet 'KEGG_RNA_DEGRADATION' with 59 genes>
    # etc.
 }
```

#### Subseting collections

Sometimes only a subset of genes is measured in an experiment. You can remove gene sets which do not contain any of the measured genes from the collection:

```python
> measured_genes = {'APOE', 'CYB5R1', 'FCER1G', 'PVR', 'HK2'}
> measured_subset = kegg_pathways.subset(measured_genes)
> print(measured_subset)
<GeneSets with 12 gene sets>
```

The skipped gene sets are accessible in `measured_subset.empty_gene_sets` for inspection.

#### Trimmming collections

```python
> kegg_pathways.trim(min_genes=10, max_genes=20)
<GeneSets with 21 gene sets>
```

#### Prettify names

```python
def prettify_kegg_name(name):
    return name.replace('KEGG_', '').replace('_', ' ')

kegg_pathways_pretty = kegg_pathways.format_names(prettify_kegg_name)
kegg_pathways_pretty.gene_sets[:2]
# (<GeneSet 'TIGHT JUNCTION' with 132 genes>, <GeneSet 'RNA DEGRADATION' with 59 genes>)
```

#### Other properties

Other properties and methods offered by `GeneSets` include:
   - `all_genes`: return a set of all genes which are covered by the gene sets in the collection
   - `name`: the name of the collection
   - `to_frame()` return a pandas `DataFrame` describing membership of the genes (gene sets = rows, genes = columns), which can be used for UpSet visualisation (e.g. with [ComplexUpset](https://github.com/krassowski/complex-upset))
   - `to_gmt(path: str)` exports the gene set to a GMT (Gene Matrix Transposed) file

## Installing GSEA implementations

Following GSEA implementations are supported:

### GSEA from Broad Institute

Login/register on [the official GSEA website](http://software.broadinstitute.org/gsea/login.jsp) and download the `gsea_3.0.jar` file (or a newer version).

Provide the location of the downloaded file to `GSEADesktop()` using `gsea_jar_path` argument, e.g.:

```python
gsea = GSEADesktop(gsea_jar_path='downloads/gsea_3.0.jar')
```

### GSEApy

To use gsea.py please install it with:

```
pip3 install gseapy
```

Use it with:

```python
from gsea_api.gsea import GSEApy

gsea = GSEApy()
```

### cudaGSEA

Please clone this fork of cudaGSEA to thirdparty directory and compile the binary version (using the instructions from [this repository](https://github.com/krassowski/cudaGSEA)):

```bash
git clone https://github.com/krassowski/cudaGSEA
```

or use [the original version](https://github.com/gravitino/cudaGSEA), which does not implement FDR calculations.

Use it with:

```python
from gsea_api.gsea import cudaGSEA

# CPU implementation can be used with use_cpu=True
gsea = cudaGSEA(fdr='full', use_cpu=False)
```

## Citation

[![DOI](https://zenodo.org/badge/188071398.svg)](https://zenodo.org/badge/latestdoi/188071398)

Please also cite the authors of the wrapped tools that you use.


## References

The initial version of this code was written for a [Master thesis project](https://github.com/krassowski/drug-disease-profile-matching) at Imperial College London.
