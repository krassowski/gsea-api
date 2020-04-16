from pytest import raises

from gsea_api.gsea.java import GSEADesktop


def test_gsea_java():
    with raises(Exception, match='Could not find GSEADesktop installation in'):
        GSEADesktop()
