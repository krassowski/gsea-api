## Legal notice

The contents of MSigDB version 6.0 and above are protected by copyright (c) 2004-2020 Broad Institute, Inc., Massachusetts Institute of Technology, and Regents of the University of California, subject to the terms and conditions of the [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).

Some MSigDB gene sets are distributed with further restrictions, see more on the [licence terms website](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).

## Test excerpts generation

The test files were generated as the tails of the original GMT files:

```bash
cat msigdb/c2.cp.reactome.v7.0.symbols.gmt | tail -n 10 > ~/gsea-api/tests/test_msigdb/c2.cp.reactome.v7.0.symbols.gmt
cat msigdb/c2.cp.reactome.v7.1.symbols.gmt | tail -n 10 > ~/gsea-api/tests/test_msigdb/c2.cp.reactome.v7.1.symbols.gmt
```
