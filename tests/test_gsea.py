from pandas import read_table, DataFrame

from gsea_api.gsea import GSEApy
from gsea_api.expression_set import ExpressionSet

matrix = read_table('tests/expression_data.tsv', index_col='Gene')
classes = ['Control'] * 4 + ['Cancer'] * 4
sub_classes = ['Control'] * 4 + ['Breast_cancer'] * 2 + ['Ovarian_cancer'] * 2


def test_gsea():
    gsea = GSEApy()
    data = ExpressionSet(matrix, classes)
    result = gsea.run(expression_data=data, **{'min-size': 1})
    assert isinstance(result, DataFrame)
    assert len(result) > 2
