from pytest import raises

from gsea_api.molecular_signatures_db import MolecularSignaturesDatabase


def test_load():
    msigdb_7_1 = MolecularSignaturesDatabase('tests/test_msigdb', version=7.1)
    assert msigdb_7_1.version == '7.1'
    assert msigdb_7_1.gene_sets == [
        {
            'name': 'c2.cp.reactome',
            'id_type': 'symbols'
        }
    ]
    reactome_7_1 = msigdb_7_1.load('c2.cp.reactome', 'symbols')
    assert 'REACTOME_NERVOUS_SYSTEM_DEVELOPMENT' in reactome_7_1.gene_sets_by_name
    assert 'REACTOME_SERINE_BIOSYNTHESIS' not in reactome_7_1.gene_sets_by_name

    msigdb_7_0 = MolecularSignaturesDatabase('tests/test_msigdb', version=7.0)
    reactome_7_0 = msigdb_7_0.load('c2.cp.reactome', 'symbols')
    assert 'REACTOME_NERVOUS_SYSTEM_DEVELOPMENT' not in reactome_7_0.gene_sets_by_name
    assert 'REACTOME_SERINE_BIOSYNTHESIS' in reactome_7_0.gene_sets_by_name


def test_fail_no_dir():
    with raises(ValueError, match='Could not find MSigDB: wrong_dir_name does not exist'):
        MolecularSignaturesDatabase('wrong_dir_name', version=7.1)
