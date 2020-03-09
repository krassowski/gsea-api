import re
from collections import Counter, defaultdict
from functools import lru_cache
from glob import glob
from pathlib import Path
from typing import Collection, Dict, List, Set
from warnings import warn

from pandas import DataFrame


class GeneSet:

    def __init__(self, name: str, genes: Collection[str], description: str = None):
        self.name = name
        self.genes = frozenset(genes)
        self.description = description

        if len(genes) == 0:
            warn(f'GeneSet {repr(name)} is empty')

        redundant_genes = None

        if len(genes) != len(self.genes):
            redundant_genes = {gene: count for gene, count in Counter(genes).items() if count > 1}

            warn(f'GeneSet {repr(name)} received a non-unique collection of genes; redundant genes: {redundant_genes}')

        self.redundant_genes = redundant_genes

    @classmethod
    def from_gmt_line(cls, line):
        name, description, *ids = line.strip().split('\t')
        return cls(name, ids, description)

    def __repr__(self):
        genes = ': ' + (', '.join(sorted(self.genes))) if len(self.genes) < 5 else ''
        return f'<GeneSet {repr(self.name)} with {len(self.genes)} genes{genes}>'


class GeneSets:

    def __init__(self, gene_sets: Collection[GeneSet], name='', allow_redundant=False):
        self.gene_sets = tuple(gene_sets)
        self.name = name
        if not allow_redundant:
            redundant = self.find_redundant()
            if any(redundant):
                message = 'Provided gene sets are not redundant; '
                if len(redundant) > 3:
                    message += (
                        f'there are {len(redundant)} gene sets having more than one name assigned; '
                        'use find_redundant() to investigate further.'
                    )
                else:
                    identical = ', '.join(
                        ' and '.join(map(repr, pathways)) + f' ({len(gene_set)} genes)'
                        for gene_set, pathways in redundant.items()
                    )
                    message += f'following gene sets are identical: {identical}'
                warn(message)

    def group_identical(self, key='name') -> Dict[frozenset, List[str]]:
        pathways_by_gene_set = defaultdict(list)
        for pathway in self.gene_sets:
            pathways_by_gene_set[pathway.genes].append(getattr(pathway, key))
        return pathways_by_gene_set

    def find_redundant(self, key='name', min_duplicates=1) -> Dict[frozenset, List[str]]:
        return {
            gene_set: pathways
            for gene_set, pathways in self.group_identical(key=key).items()
            if len(pathways) > min_duplicates
        }

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
    @lru_cache()
    def all_genes(self):
        all_genes = set()

        for gene_set in self.gene_sets:
            all_genes.update(gene_set.genes)

        return all_genes

    def to_frame(self) -> DataFrame:
        all_genes = self.all_genes
        return DataFrame(
            [
                [
                    gene in gene_set.genes
                    for gene in all_genes
                ]
                for gene_set in self.gene_sets
            ],
            index=[gene_set.name for gene_set in self.gene_sets],
            columns=all_genes
        )

    @property
    @lru_cache()
    def gene_sets_by_name(self):
        by_name = {
            gene_set.name: gene_set
            for gene_set in self.gene_sets
        }
        assert len(self.gene_sets) == len(by_name)
        return by_name

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
