from pandas import read_table, DataFrame
from pytest import raises

from gsea_api import cudaGSEA
from gsea_api.gsea import GSEApy
from gsea_api.expression_set import ExpressionSet
from gsea_api.gsea.exceptions import GSEANoResults
from gsea_api.molecular_signatures_db import GeneSets

matrix = read_table('tests/expression_data.tsv', index_col='Gene')
classes = ['Control'] * 4 + ['Cancer'] * 4


def test_gsea_py():
    gsea = GSEApy()
    data = ExpressionSet(matrix.copy(), classes)
    result = gsea.run(expression_data=data, **{'min-size': 1})
    assert isinstance(result, DataFrame)
    assert len(result) > 2


def test_cuda_gsea(monkeypatch):
    gsea = cudaGSEA(path='cudaGSEA/cudaGSEA/src/cudaGSEA')
    data = ExpressionSet(matrix.copy(), classes)
    gene_sets = GeneSets.from_gmt('tests/sets_for_expression.gmt')

    def mock_cuda_run(*args, **kwargs):
        return read_table('tests/cudaGSEA_subprocess_output.tsv')

    monkeypatch.setattr(cudaGSEA, 'run_in_subprocess', mock_cuda_run)

    result = gsea.run(expression_data=data, gene_sets=gene_sets)
    assert isinstance(result, DataFrame)
    assert len(result) > 2

    def mock_cuda_run_no_results(*args, **kwargs):
        return DataFrame()

    monkeypatch.setattr(cudaGSEA, 'run_in_subprocess', mock_cuda_run_no_results)
    data = ExpressionSet(matrix.copy(), classes)

    with raises(GSEANoResults):
        gsea.run(expression_data=data, gene_sets=gene_sets)
