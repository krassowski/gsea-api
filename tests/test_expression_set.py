from pandas import read_table, Series
from pytest import raises

from gsea_api.expression_set import ExpressionSet

matrix = read_table('tests/expression_data.tsv', index_col='Gene')
classes = ['Control'] * 4 + ['Cancer'] * 4
sub_classes = ['Control'] * 4 + ['Breast_cancer'] * 2 + ['Ovarian_cancer'] * 2


def test_from_frame():

    # accepts list classes:
    data = ExpressionSet(matrix, classes)
    assert all(data.classes == Series(classes))

    # accepts Series classes:
    data = ExpressionSet(matrix, Series(classes))
    assert all(data.classes == Series(classes))

    with raises(ValueError, match='Number of classes different from the number of columns'):
        first_five_patients = matrix.iloc[:, :5]
        ExpressionSet(first_five_patients, classes)


def test_contrast():
    data = ExpressionSet(matrix, sub_classes)
    subset = data.contrast(case='Breast_cancer', control='Control')
    assert isinstance(subset, ExpressionSet)
    assert len(subset.joined) == len(data.joined)
    assert len(subset.classes) != len(data.classes)
    assert len(subset.classes) == 4 + 2
    assert set(subset.classes) == {'Breast_cancer', 'Control'}


def test_hashable():
    data = ExpressionSet(matrix, sub_classes)
    subset = data.contrast(case='Breast_cancer', control='Control')
    assert hash(data.hashable) != hash(subset.hashable)
