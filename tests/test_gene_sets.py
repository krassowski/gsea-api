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
    # max_genes should be optional
    assert len(pathways.trim(min_genes=3)) == 2
    # min_genes should be optional
    assert len(pathways.trim(max_genes=2)) == 2

    assert pathways.name == 'gene_ontology_pathways.gmt'
    assert repr(pathways) == "<GeneSets 'gene_ontology_pathways.gmt' with 4 gene sets>"

    uterus_morphogenesis = pathways.gene_sets_by_name['GO:0061038']
    assert uterus_morphogenesis.description == 'uterus morphogenesis'
    assert uterus_morphogenesis.genes == frozenset('ASH1L	KDM5B	STRA6	WNT7A	WNT9B	NIPBL'.split())


def test_to_frame():
    gene_sets = GeneSets(gene_sets_list)
    df = gene_sets.to_frame()
    assert len(df) == 3
    assert set(df.columns) == set(all_genes)
    assert df[all_genes].equals(
        read_table(StringIO(expected_frame), index_col='gene_set').astype(bool)
    )


def test_gene_sets():
    identical_gene_sets = r"'behavioral response to chemical pain' and 'behavioral response to formalin induced pain' \(2 genes\)"

    warning_text = f'Provided gene sets are not redundant; following gene sets are identical: {identical_gene_sets}'
    with warns(UserWarning, match=warning_text):
        gene_sets = GeneSets(gene_sets_list)
    assert len(gene_sets) == len(gene_sets_list)
    assert len(gene_sets) == 3

    assert gene_sets.all_genes == set(all_genes)
    assert repr(gene_sets) == "<GeneSets with 3 gene sets>"

    subset = gene_sets.subset({'ASHL1', 'KDM5B'})
    assert subset != gene_sets
    assert list(subset.gene_sets_by_name.keys()) == ['uterus morphogenesis']

    # subset() should accept any iterable, not only sets
    assert subset == gene_sets.subset(['ASHL1', 'KDM5B'])

    warning_text = f'Provided gene sets are not redundant; there are 10 gene sets having more than one name assigned'
    with warns(UserWarning, match=warning_text):
        gene_sets = GeneSets([
            GeneSet('duplicated set', [f'G{i}', f'H{i}'])
            for i in range(10)
            for _ in range(2)
        ])

    with warns(UserWarning, match='There are 2 empty gene sets.*'):
        gene_sets = GeneSets([
            GeneSet('empty gene set 1', genes=[], warn_if_empty=False),
            GeneSet('empty gene set 2', genes=[], warn_if_empty=False)
        ])
        assert len(gene_sets) == 0
        assert len(gene_sets.empty_gene_sets) == 2


def test_gene_set():
    gene_set = GeneSet('behavioral response to chemical pain', ['P2RX3', 'NTRK1'])
    assert repr(gene_set) == "<GeneSet 'behavioral response to chemical pain' with 2 genes: NTRK1, P2RX3>"

    with warns(UserWarning, match="GeneSet 'non-unique collection' received a non-unique collection of genes; redundant genes: {'TP53': 2}"):
        gene_set = GeneSet('non-unique collection', ['TP53', 'TP53'])
    assert gene_set.genes == {'TP53'}

    with warns(UserWarning, match="GeneSet 'empty collection' is empty"):
        gene_set = GeneSet('empty collection', [])

    with warns(None) as record:
        gene_set = GeneSet('empty collection', [], warn_if_empty=False)
    assert not record.list

