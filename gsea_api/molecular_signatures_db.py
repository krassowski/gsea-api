import re
from glob import glob
from pathlib import Path
from typing import Set, Collection
from warnings import warn

from pandas import DataFrame


class GeneSet:

    def __init__(self, name: str, genes: Collection[str]):
        self.name = name
        self.genes = set(genes)

    @classmethod
    def from_gmt_line(cls, line):
        name, url, *ids = line.strip().split('\t')
        return cls(name, ids)

    def __repr__(self):
        genes = ': ' + (', '.join(sorted(self.genes))) if len(self.genes) < 5 else ''
        return f'<GeneSet {repr(self.name)} with {len(self.genes)} genes{genes}>'


class GeneSets:

    def __init__(self, gene_sets: Collection[GeneSet], name='', allow_redundant=False):
        self.gene_sets = gene_sets
        self.name = name
        if not allow_redundant:
            df = self.to_frame()
            if any(df.duplicated()):
                identical = ' AND '.join(df.index[df.duplicated(keep=False)])
                warn(f'Provided gene sets are not redundant; following gene sets are identical: {identical}')

    @classmethod
    def from_gmt(cls, path, name=None):
        with open(path) as f:
            return cls({
                GeneSet.from_gmt_line(line)
                for line in f
            }, name=name or Path(path).name)

    def trim(self, min_genes, max_genes: int):
        return GeneSets({
            gene_set
            for gene_set in self.gene_sets
            if min_genes <= len(gene_set.genes) <= max_genes
        })

    def format_names(self, formatter):
        for gene_set in self.gene_sets:
            gene_set.name = formatter(gene_set.name)
        return self

    def _to_gmt(self, f):
        for gene_set in self.gene_sets:
            f.write(gene_set.name + '\t' + '\t'.join(gene_set.genes) + '\n')

    def to_gmt(self, path):
        if isinstance(path, str):
            with open(path, mode='w') as f:
                self._to_gmt(f)
        else:
            self._to_gmt(path)

    def subset(self, genes: Set[str]):
        return GeneSets({
            GeneSet(name=gene_set.name, genes=gene_set.genes & genes)
            for gene_set in self.gene_sets
        })

    @property
    def all_identifiers(self):
        all_identifiers = set()

        for gene_set in self.gene_sets:
            all_identifiers.update(gene_set.genes)

        return all_identifiers

    def to_frame(self) -> DataFrame:
        identifiers = self.all_identifiers
        return DataFrame(
            [
                [
                    gene in gene_set.genes
                    for gene in identifiers
                ]
                for gene_set in self.gene_sets
            ],
            index=[gene_set.name for gene_set in self.gene_sets],
            columns=identifiers
        )

    def __len__(self):
        return len(self.gene_sets)

    def __repr__(self):
        name = ' ' + repr(self.name) if self.name else ''
        return f'<GeneSets{name} with {len(self.gene_sets)} gene sets>'


# for backwards compatibility
GeneMatrixTransposed = GeneSets


class MolecularSignaturesDatabase:
    def __init__(self, path, version='6.2'):
        self.path = Path(path)
        self.version = version
        wildcard_path = str((self.path / f'*.v{self.version}.*.gmt').resolve())
        self.gene_sets = [
            self.parse_name(Path(p).name)
            for p in glob(wildcard_path)
        ]

    def parse_name(self, name):
        parsed = re.match(rf'(?P<name>.*?)\.v{self.version}\.(?P<id_type>(entrez|symbols)).gmt', name)
        return parsed.groupdict()

    def resolve(self, gene_sets, id_type):
        path = self.path / f'{gene_sets}.v6.2.{id_type}.gmt'
        if path.exists():
            return str(path)
        else:
            raise ValueError(f'Unknown library: {path}!')

    def load(self, gene_sets, id_type) -> GeneSets:
        path = self.resolve(gene_sets=gene_sets, id_type=id_type)

        return GeneSets.from_gmt(path, name=gene_sets)
