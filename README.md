# GSEA API for Pandas
Pandas API for Gene Set Enrichment Analysis in Python (GSEApy, cudaGSEA, GSEA)

- This Python wrapper around various GSEA implementations aims to provide a unified programming interface,
built using the pandas DataFrames and a hierarchy of Pythonic classes.
- The file exports (providing input for GSEA) were written with performance in mind, using lower level numpy functions where necessary, thus are much faster than usual pandas-based exports.
- This project aims to allow scientists in the Python community to easily compare different implementations of GSEA, and to integrate those in projects which require high performance GSEA interface.
- The project is in work-in-progress state and scheduled to have a major refactor and a more complete documentation.

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

#### Installing cudaGSEA

Please clone this fork of cudaGSEA to thirdparty directory and compile the binary version:

```
git clone https://github.com/krassowski/cudaGSEA
```

or use [the original version](https://github.com/gravitino/cudaGSEA), which does not implement FDR calculations.

### Citation

Please cite the authors of the wrapped tools that you use.


### References

The initial version of this code was written for my [Master thesis project](https://github.com/krassowski/drug-disease-profile-matching) at Imperial College London.

