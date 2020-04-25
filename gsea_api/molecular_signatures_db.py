import re
from collections import Counter, defaultdict
from functools import lru_cache
from glob import glob
from pathlib import Path
from typing import Collection, Dict, List, Iterable
from warnings import warn

from pandas import DataFrame


class GeneSet:

    def __init__(self, name: str, genes: Collection[str], description: str = None, warn_if_empty=True):
        self.name = name
        self.genes = frozenset(genes)
        self.description = description

        if warn_if_empty and self.is_empty:
            warn(f'GeneSet {repr(name)} is empty')

        redundant_genes = None

        if len(genes) != len(self.genes):
            redundant_genes = {gene: count for gene, count in Counter(genes).items() if count > 1}

            warn(f'GeneSet {repr(name)} received a non-unique collection of genes; redundant genes: {redundant_genes}')

        self.redundant_genes = redundant_genes

    @classmethod
    def from_gmt_line(cls, line, **kwargs):
        name, description, *ids = line.strip().split('\t')
        return cls(name, ids, description, **kwargs)

    @property
    def is_empty(self):
        return len(self.genes) == 0

    def __repr__(self):
        genes = ': ' + (', '.join(sorted(self.genes))) if len(self.genes) < 5 else ''
        return f'<GeneSet {repr(self.name)} with {len(self.genes)} genes{genes}>'

    def __eq__(self, other: 'GeneSet'):
        return (
            self.name == other.name
            and
            self.genes == other.genes
        )

    def __hash__(self):
        return hash((self.name, self.genes))


class GeneSets:

    def __init__(self, gene_sets: Collection[GeneSet], name='', allow_redundant=False, remove_empty=True):
        self.gene_sets = tuple(gene_sets)   # NOTE: this is not final
        self.name = name
        if not allow_redundant:
            redundant = self.find_redundant()
            if any(redundant):
                message = 'Provided gene sets are not redundant; '
                if len(redundant) > 3:
                    message += (
                        f'there are {len(redundant)} gene sets having more than one name assigned; '
                        'use `find_redundant()` to investigate further.'
                    )
                else:
                    identical = ', '.join(
                        ' and '.join(map(repr, pathways)) + f' ({len(gene_set)} genes)'
                        for gene_set, pathways in redundant.items()
                    )
                    message += f'following gene sets are identical: {identical}'
                warn(message)

        empty_gene_sets = {gene_set for gene_set in gene_sets if gene_set.is_empty}

        if len(empty_gene_sets) != 0:
            empty_message = (
                ', '.join(gene_set.name for gene_set in empty_gene_sets)
                if len(empty_gene_sets) <= 5 else
                'use `empty_gene_sets` property to investigate further.'
            )
            warn(f'There are {len(empty_gene_sets)} empty gene sets: {empty_message}')

            if remove_empty:
                gene_sets = set(gene_sets) - empty_gene_sets
                warn(f'{len(empty_gene_sets)} empty gene sets were removed.')

        self.empty_gene_sets = empty_gene_sets
        self.gene_sets = tuple(gene_sets)

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
    def from_gmt(cls, path, name=None, **kwargs):
        with open(path) as f:
            return cls(
                {
                    GeneSet.from_gmt_line(line, warn_if_empty=False)
                    for line in f
                },
                name=name or Path(path).name,
                **kwargs
            )

    def trim(self, min_genes: int = 0, max_genes: int = float('Inf')):
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

    def subset(self, genes: Iterable[str]):
        if not isinstance(genes, set):
            genes = set(genes)
        return GeneSets({
            GeneSet(name=gene_set.name, genes=gene_set.genes & genes, warn_if_empty=False)
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

    def __eq__(self, other: 'GeneSets'):
        return (
            set(self.gene_sets) == set(other.gene_sets)
            and
            self.name == other.name
        )

    def __hash__(self):
        return hash((self.name, self.gene_sets))


# for backwards compatibility
GeneMatrixTransposed = GeneSets


class MolecularSignaturesDatabase:
    def __init__(self, path, version='7.1'):
        self.path = Path(path)
        if not self.path.exists():
            raise ValueError(f'Could not find MSigDB: {self.path} does not exist')
        self.version = str(version)
        wildcard_path = str((self.path / f'*.v{self.version}.*.gmt').resolve())
        self.gene_sets = [
            self.parse_name(Path(p).name)
            for p in glob(wildcard_path)
        ]

    def parse_name(self, name):
        parsed = re.match(rf'(?P<name>.*?)\.v{self.version}\.(?P<id_type>(entrez|symbols)).gmt', name)
        return parsed.groupdict()

    def resolve(self, gene_sets, id_type):
        path = self.path / f'{gene_sets}.v{self.version}.{id_type}.gmt'
        if path.exists():
            return str(path)
        else:
            raise ValueError(f'Unknown library: {path}!')

    def load(self, gene_sets, id_type) -> GeneSets:
        path = self.resolve(gene_sets=gene_sets, id_type=id_type)

        return GeneSets.from_gmt(path, name=gene_sets)
