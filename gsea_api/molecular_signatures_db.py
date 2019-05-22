import re
from glob import glob
from pathlib import Path
from typing import Set


class GeneSet:

    def __init__(self, name, genes):
        self.name = name
        self.genes = set(genes)

    @classmethod
    def from_gmt_line(cls, line):
        name, url, *ids = line.strip().split('\t')
        return cls(name, ids)


class GeneMatrixTransposed:

    def __init__(self, gene_sets, name=''):
        self.gene_sets = gene_sets
        self.name = name

    @classmethod
    def from_gmt(cls, path, name=''):
        with open(path) as f:
            return cls({
                GeneSet.from_gmt_line(line)
                for line in f
            }, name=name)

    def trim(self, min_genes, max_genes: int):
        return GeneMatrixTransposed({
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
        return GeneMatrixTransposed({
            GeneSet(name=gene_set.name, genes=gene_set.genes & genes)
            for gene_set in self.gene_sets
        })

    @property
    def all_identifiers(self):
        all_identifiers = set()

        for gene_set in self.gene_sets:
            all_identifiers.update(gene_set.genes)

        return all_identifiers


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

    def load(self, gene_sets, id_type) -> GeneMatrixTransposed:
        path = self.resolve(gene_sets=gene_sets, id_type=id_type)

        return GeneMatrixTransposed.from_gmt(path, name=gene_sets)
