import pandas
from numpy import savetxt
from pandas import DataFrame, Series, concat
from typing.io import TextIO
from typing import List, Union


class ExpressionSet:
    # TODO: have a limit() method (to trim the data to the most relevant genes only)

    joined: DataFrame
    classes: Series

    def __init__(self, data: DataFrame, classes: Union[Series, List]):
        self.joined = data
        self.classes = Series(classes)
        if len(data.columns) != len(classes):
            raise ValueError('Number of classes different from the number of columns')

    @classmethod
    def from_cases_and_controls(
        cls, cases: DataFrame, controls: DataFrame,
        case_name: str = 'case', control_name: str = 'control'
    ) -> 'ExpressionSet':
        assert cases.index == controls.index
        return cls(
            data=concat(cases, controls),
            classes=Series(
                [case_name] * len(cases.columns) + [control_name] * len(controls.columns)
            )
        )

    #def __repr__(self):
    #    return f'<{self.__class__.__name__}: {len(self.cases.columns)} cases, {len(self.controls.columns)} controls>'

    def contrast(self, case, control) -> 'ExpressionSet':
        """Select a subset of two groups (cases and controls) from experimental groups in the expression set.

        For example, having groups Drug_A, Drug_B, Control, if you want to only analyse Drug_A vs Control, specify:

        subset = expression_set.contrast(case='Drug_A', control='Control')
        """
        assert case in self.classes.unique() and control in self.classes.unique()
        mask = self.classes.isin({case, control})
        return ExpressionSet(
            self.joined[self.joined.columns[mask]],
            self.classes[mask]
        )

    @property
    def hashable(self):
        all_data = self.joined
        return tuple(all_data.index), tuple(all_data.columns), all_data.sum().sum()

    @property
    def safe_classes(self):
        return [cls.replace(' ', '_') for cls in self.classes]

    def to_cls(self, f):
        classes = self.safe_classes
        f.write(f'{len(classes)} {len(set(classes))} 1\n')
        classes_set = []
        for class_ in classes:
            if class_ not in classes_set:
                classes_set.append(class_)
        f.write(f'# {" ".join(classes_set)}\n')
        f.write(' '.join(classes))

    def to_gct(self, f: TextIO, tabular_writer='to_txt'):
        f.write('#1.2\n')
        expression_data = self.joined
        assert expression_data.notnull().all().all()
        f.write(f'{len(expression_data)}\t{len(expression_data.columns)}\n')
        getattr(self, tabular_writer)(f, expression_data)

    def to_txt(self, f: TextIO, expression_data=None):
        if expression_data is None:
            from copy import copy
            expression_data = copy(self.joined)
        columns = expression_data.columns
        pandas.options.mode.chained_assignment = None
        expression_data['Description'] = 'na'
        expression_data = expression_data[['Description', *columns]]
        if type(expression_data.index[0]) is bytes:
            expression_data.index = [b.decode('utf-8') for b in expression_data.index]
        expression_data.index = expression_data.index.astype(str)
        expression_data.index.name = 'gene'
        header = '\t'.join([expression_data.index.name, *expression_data.columns]) + '\n'
        f.write(header)

        savetxt(
            f,
            expression_data.reset_index().values,
            delimiter='\t',
            # entrez id or symbol (as str), 'na' (not a description, str), *data (floats)
            fmt='%s\t%s' + ('\t%f' * (len(expression_data.columns) - 1))
        )
