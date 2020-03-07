# GSEA API for Pandas
Pandas API for Gene Set Enrichment Analysis in Python (GSEApy, cudaGSEA, GSEA)

- This Python wrapper aims to provide a unified API for various GSEA implementations; it uses pandas DataFrames and a hierarchy of Pythonic classes.
- The file exports (providing input for GSEA) were written with performance in mind, using lower level numpy functions where necessary, thus are much faster than pandas-based exports.
- This project aims to allow researchers to easily compare different implementations of GSEA, and to integrate those in projects which require high performance GSEA.
- The project is in work-in-progress state and may undergo moderate refactoring and recieve a more complete documentation (if there is an interest in it).

### Example usage

```python
from pandas import read_csv
from gsea_api.expression_set import ExpressionSet
from gsea_api.gsea import GSEADesktop
from gsea_api.molecular_signatures_db import GeneMatrixTransposed

reactome_pathways = GeneMatrixTransposed.from_gmt('ReactomePathways.gmt')

gsea = GSEADesktop()

design = ['Disease', 'Disease', 'Disease', 'Control', 'Control', 'Control']
matrix = read_csv('expression_data.csv')

result = gsea.run(
    # note: contrast() is not necessary in this simple case
    ExpressionSet(matrix, design).contrast('Disease', 'Control'),
    reactome_pathways,
    metric='Signal2Noise',
    permutations=1000
)
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

The initial version of this code was written for my [Master thesis project](https://github.com/krassowski/drug-disease-profile-matching) at Imperial College London.

