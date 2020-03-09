from io import StringIO
from pandas import read_table

from pytest import warns
from gsea_api.molecular_signatures_db import GeneSets, GeneSet

gene_sets_list = [
    GeneSet('uterus morphogenesis', 'ASH1L	KDM5B	STRA6	WNT7A	WNT9B	NIPBL'.split()),
    GeneSet('behavioral response to chemical pain', ['P2RX3', 'NTRK1']),
    GeneSet('behavioral response to formalin induced pain', ['P2RX3', 'NTRK1'])
]
all_genes = 'ASH1L	KDM5B	STRA6	WNT7A	WNT9B	NIPBL	P2RX3	NTRK1'.split()


expected_frame = """
gene_set	ASH1L	KDM5B	STRA6	WNT7A	WNT9B	NIPBL	P2RX3	NTRK1
uterus morphogenesis	1	1	1	1	1	1	0	0
behavioral response to chemical pain	0	0	0	0	0	0	1	1
behavioral response to formalin induced pain	0	0	0	0	0	0	1	1
"""


def test_from_gmt():
    pathways = GeneSets.from_gmt('tests/gene_ontology_pathways.gmt')
    assert len(pathways) == len(pathways.gene_sets)
    assert len(pathways.gene_sets) == 4

    assert len(pathways.trim(min_genes=3, max_genes=100)) == 2
    assert pathways.name == 'gene_ontology_pathways.gmt'
    assert repr(pathways) == "<GeneSets 'gene_ontology_pathways.gmt' with 4 gene sets>"

    uterus_morphogenesis = pathways.gene_sets_by_name['GO:0061038']
    assert uterus_morphogenesis.description == 'uterus morphogenesis'
    assert uterus_morphogenesis.genes == set('ASH1L	KDM5B	STRA6	WNT7A	WNT9B	NIPBL'.split())


def test_to_frame():
    gene_sets = GeneSets(gene_sets_list)
    df = gene_sets.to_frame()
    assert len(df) == 3
    assert set(df.columns) == set(all_genes)
    assert df[all_genes].equals(
        read_table(StringIO(expected_frame), index_col='gene_set').astype(bool)
    )


def test_gene_sets():
    identical_gene_sets = ' AND '.join([
        'behavioral response to chemical pain',
        'behavioral response to formalin induced pain'
    ])
    warning_text = f'Provided gene sets are not redundant; following gene sets are identical: {identical_gene_sets}'
    with warns(UserWarning, match=warning_text):
        gene_sets = GeneSets(gene_sets_list)
    assert len(gene_sets) == len(gene_sets_list)
    assert len(gene_sets) == 3

    assert gene_sets.all_genes == set(all_genes)
    assert repr(gene_sets) == "<GeneSets with 3 gene sets>"


def test_gene_set():
    gene_set = GeneSet('behavioral response to chemical pain', ['P2RX3', 'NTRK1'])
    assert repr(gene_set) == "<GeneSet 'behavioral response to chemical pain' with 2 genes: NTRK1, P2RX3>"
