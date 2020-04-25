from copy import copy
from os import makedirs
from os import remove
from pathlib import Path
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from typing import Set, Union
from warnings import warn

from pandas import DataFrame
from jupyter_helpers.following_tail import FollowingTail

from abc import ABC, abstractmethod, abstractproperty

from ..molecular_signatures_db import GeneSets
from ..streams import handle_streams
from ..expression_set import ExpressionSet

from .exceptions import GSEANotEnoughSamples
from ..paths import tmp_dir


class GSEA(ABC):

    def __init__(self, display=None, display_error=warn):
        self._original_display = display
        self.display = self.get_display(display)
        self.display_error = display_error

    @staticmethod
    def get_display(display):
        if not display:
            display = FollowingTail()
        return display

    def prepare_output(self):
        self.display = self.get_display(self._original_display)
        if isinstance(self.display, FollowingTail):
            self.display.activate()

    @staticmethod
    def prepare_files(expression_data: ExpressionSet):
        with NamedTemporaryFile('w', delete=False, dir=tmp_dir, suffix='.cls') as f:
            classes_path = f.name
            expression_data.to_cls(f)

        with NamedTemporaryFile('w', delete=False, dir=tmp_dir, suffix='.gct') as f:
            expression_data.to_gct(f)
            data_path = f.name

        return data_path, classes_path

    @property
    def default_handlers(self):
        return {
            'out': self.display,
            'err': self.display_error
        }

    def read_from_process(self, process, handlers=None, verbose=False):
        if handlers is None:
            handlers = copy(self.default_handlers)

        if not verbose:
            del self.default_handlers['out']

        return handle_streams(process, handlers)

    def run_in_subprocess(self, command, verbose=False):
        if verbose:
            print(command)
        arguments = command.split(' ')
        process = Popen(arguments, stdout=PIPE, stderr=PIPE)
        return self.read_from_process(process, verbose=verbose)

    @staticmethod
    def prepare_out_dir(out_dir, temp_dir):
        if not out_dir:
            out_dir = Path(temp_dir)
        else:
            out_dir = Path(out_dir)
            if not out_dir.is_absolute():
                pwd = Path().absolute().parent
                out_dir = pwd / str(out_dir)
            makedirs(out_dir, exist_ok=True)
        return out_dir

    @abstractproperty
    def metrics_using_variance(self) -> Set[str]:
        """Metrics which would fail if given less than three samples per category"""

    def has_enough_samples_for_metric(self, expression_data, metric: str):
        if metric in self.metrics_using_variance:
            # TODO
            return True
            #return len(expression_data.controls.columns) >= 3 and len(expression_data.cases.columns) >= 3
        return True

    @abstractmethod
    def run(
        self,
        expression_data: ExpressionSet, gene_sets: Union[str, GeneSets],
        metric: str = None, **kwargs
    ) -> DataFrame:
        """Run the GSEA method using case and control data from given expression object.

        Raises:
            GSEANoResults exception if the method fails with no results.
            GSEANotEnoughSamples if specified metric requires more samples than provided.
        """

        if not self.has_enough_samples_for_metric(expression_data, metric):
            warn(f'Too few samples for the metric {metric}')
            raise GSEANotEnoughSamples(f'Too few samples for the metric {metric}')

    @staticmethod
    def clean_up(data_path, classes_path):
        remove(data_path)
        remove(classes_path)

    @staticmethod
    def solve_gene_sets(gene_sets):
        if isinstance(gene_sets, GeneSets):
            if gene_sets.path:
                gene_sets = gene_sets.path
            else:
                with NamedTemporaryFile(mode='w', delete=False, suffix= '.gmt') as f:
                    gene_sets.to_gmt(f)
                    gene_sets = f.name
        return gene_sets
