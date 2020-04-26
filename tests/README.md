## Legal

A subset of Gene Ontology gene sets is used in tests. Gene Ontology is licenced under [Creative Commons Attribution 4.0 Unported License](https://creativecommons.org/licenses/by/4.0/legalcode).

```python
from gsea_api.molecular_signatures_db import MolecularSignaturesDatabase
(
    MolecularSignaturesDatabase('msigdb', version=7.1)
    .load('c5.bp', 'symbols')
    .subset(['TACC2', 'TP53', 'BRCA1', 'NOTCH1'])
    .trim(3)
    .to_gmt('sets_for_expression.gmt')
)
```
