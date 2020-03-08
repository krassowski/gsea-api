# GSEA API for Pandas
[![Build Status](https://travis-ci.com/krassowski/gsea-api.svg?branch=master)](https://travis-ci.com/krassowski/gsea-api) 
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](http://choosealicense.com/licenses/mit/)
[![DOI](https://zenodo.org/badge/188071398.svg)](https://zenodo.org/badge/latestdoi/188071398)

Pandas API for Gene Set Enrichment Analysis in Python (GSEApy, cudaGSEA, GSEA)

- aims to provide a unified API for various GSEA implementations; uses pandas DataFrames and a hierarchy of Pythonic classes.
- file exports (exporting input for GSEA) use low-level numpy functions and are much faster than in pandas
- aims to allow researchers to easily compare different implementations of GSEA, and to integrate those in projects which require high-performance GSEA (e.g. massive screening for drug-repositioning)
- provides useful utilities for work with GMT files, or gene sets and pathways in general in Python

### Example usage

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

### Installation

To install the API use:
```
pip3 install gsea_api
```

#### Installing GSEA from Broad Institute

Login/register on [the official GSEA website](http://software.broadinstitute.org/gsea/login.jsp) and download the `gsea_3.0.jar` file (or a newer version).

Please place the downloaded file in the thirdparty directory.


#### Installing GSEApy

To use gsea.py please install it with:

```
pip3 install gseapy
```

and link its binary to the `thirdparty` directory
```
ln -s virtual_environment_path/bin/gseapy thirdparty/gseapy
```


Use it with:

```python
from gsea_api.gsea import GSEApy

gsea = GSEApy()
```

#### Installing cudaGSEA

Please clone this fork of cudaGSEA to thirdparty directory and compile the binary version (using the instructions from [this repository](https://github.com/krassowski/cudaGSEA)):

```
git clone https://github.com/krassowski/cudaGSEA
```

or use [the original version](https://github.com/gravitino/cudaGSEA), which does not implement FDR calculations.

Use it with:

```python
from gsea_api.gsea import cudaGSEA

# CPU implementation can be used with use_cpu=True
gsea = cudaGSEA(fdr='full', use_cpu=False)
```

### Citation

[![DOI](https://zenodo.org/badge/188071398.svg)](https://zenodo.org/badge/latestdoi/188071398)

Please also cite the authors of the wrapped tools that you use.


### References

The initial version of this code was written for a [Master thesis project](https://github.com/krassowski/drug-disease-profile-matching) at Imperial College London.
