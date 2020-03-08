from gsea_api.molecular_signatures_db import GeneSets


def test_from_gmt():
    pathways = GeneSets.from_gmt('tests/gene_ontology_pathways.gmt')
    assert len(pathways) == len(pathways.gene_sets)
    assert len(pathways.gene_sets) == 4

    assert len(pathways.trim(min_genes=3, max_genes=100)) == 2
