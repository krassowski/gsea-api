import re
from collections import Counter, defaultdict
from copy import deepcopy
from functools import lru_cache
from glob import glob
from pathlib import Path
from typing import Collection, Dict, List, Iterable, Callable, Union, TextIO, Set, TypeVar
from warnings import warn
from xml.etree import ElementTree

from pandas import DataFrame


class GeneSet:

    def __init__(
        self,
        name: str,
        genes: Collection[str],
        description: str = None,
        warn_if_empty: bool = True,
        representativeness: float = None,
        metadata: Dict[str, str] = None
    ):
        self.name = name
        self.genes = frozenset(genes)
        self.description = description
        self.representativeness = representativeness
        self.metadata = metadata or {}

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
    def is_empty(self) -> bool:
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


T = TypeVar('T')


def lru_cache_mypy(func: Callable[..., T]) -> T:
    """Workaround for mypy issue with lru_cache on property see:
    https://github.com/python/mypy/issues/5858#issuecomment-454144705
    """
    return lru_cache()(func)  # type: ignore


class GeneSets:

    def __init__(
        self,
        gene_sets: Collection[GeneSet],
        name: str = '',
        allow_redundant: bool = False,
        remove_empty: bool = True,
        path: Union[str, Path] = None,
        collapse_redundant: Union[str, bool] = False,
        collapse_limit: int = 10
    ):
        self.gene_sets = tuple(gene_sets)   # NOTE: this is not final
        self.name = name
        self.path = path

        redundant = self.find_redundant()
        redundant_gene_sets_n = len(sum(redundant.values(), []))
        self.redundant = redundant

        if not allow_redundant:
            # do not warn if we will be collapsing anyways
            if any(redundant) and not collapse_redundant:
                message = 'Provided gene sets are not redundant; '
                if redundant_gene_sets_n > 3:
                    message += (
                        f'there are {len(redundant)} gene sets having more than one name assigned,'
                        f' in total affecting {redundant_gene_sets_n} gene sets'
                    )
                else:
                    identical = ', '.join(
                        ' and '.join(map(repr, pathways)) + f' ({len(gene_set)} genes)'
                        for gene_set, pathways in redundant.items()
                    )
                    message += f'following gene sets are identical: {identical}'
                message += (
                    '; use `.find_redundant()` to investigate further,'
                    ' or collapse them by passing `collapse_redundant: str`.'
                )
                warn(message)

        empty_gene_sets = {gene_set for gene_set in gene_sets if gene_set.is_empty}

        if len(empty_gene_sets) != 0:
            if remove_empty:
                gene_sets = set(gene_sets) - empty_gene_sets
                warn(f'{len(empty_gene_sets)} empty gene sets were removed.')
            else:
                empty_message = (
                    ', '.join(gene_set.name for gene_set in empty_gene_sets)
                    if len(empty_gene_sets) <= 5 else
                    'use `empty_gene_sets` property to investigate further.'
                )
                warn(f'There are {len(empty_gene_sets)} empty gene sets: {empty_message}')

        if collapse_redundant and any(redundant):
            redundant_gene_frozen_sets = {
                frozen_genes
                for frozen_genes in redundant.keys()
            }
            non_redundant_gene_sets = [
                gene_set
                for gene_set in gene_sets
                if frozenset(gene_set.genes) not in redundant_gene_frozen_sets
            ]
            for frozen_genes, names in redundant.items():
                if len(names) > collapse_limit:
                    examples = ','.join(names[:3])
                    warn(
                        f'Redundant collection of {len(names)} gene sets which includes {examples} (among others)'
                        f' exceeds the provided `collapse_limit` ({collapse_limit}) and only as many'
                        f' names will be included in the collapse gene set.'
                    )
            collapsed_redundant_gene_sets = [
                GeneSet(
                    genes=frozen_genes,
                    name=collapse_redundant.join(names[:collapse_limit]) + (
                        f'{collapse_redundant}... {len(names) - collapse_limit} more'
                        if len(names) > collapse_limit else
                        ''
                    ),
                    description=collapse_redundant.join(names)
                )
                for frozen_genes, names in redundant.items()
            ]
            warn(
                f'Collapsed {redundant_gene_sets_n} redundant gene sets'
                f' into {len(redundant)} non-redundant sets,'
                f' concatenating their names with {repr(collapse_redundant)}'
            )
            gene_sets = collapsed_redundant_gene_sets + non_redundant_gene_sets

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
                path=path,
                **kwargs
            )

    def trim(self, min_genes: int = 0, max_genes: int = float('Inf')) -> 'GeneSets':
        return GeneSets({
            gene_set
            for gene_set in self.gene_sets
            if min_genes <= len(gene_set.genes) <= max_genes
        })

    def collapse_redundant(self, sep: str) -> 'GeneSets':
        return GeneSets(self.gene_sets, collapse_redundant=sep)

    def format_names(self, formatter: Callable[[GeneSet], str]) -> 'GeneSets':
        result = deepcopy(self)
        for gene_set in result.gene_sets:
            gene_set.name = formatter(gene_set)
        return result

    def _to_gmt(self, f: TextIO):
        for gene_set in self.gene_sets:
            f.write(gene_set.name + '\t' + '\t'.join(gene_set.genes) + '\n')

    def to_gmt(self, path: Union[TextIO, str, Path]):
        if isinstance(path, str) or isinstance(path, Path):
            with open(path, mode='w') as f:
                self._to_gmt(f)
        else:
            self._to_gmt(path)

    def extract(self, set_names: Iterable[str]) -> 'GeneSets':
        set_names = set(set_names)
        return GeneSets({
            gene_set
            for gene_set in self.gene_sets
            if gene_set.name in set_names
        })

    def subset(self, genes: Iterable[str], min_representation: float = 0) -> 'GeneSets':
        if not isinstance(genes, set):
            genes = set(genes)
        return GeneSets({
            GeneSet(
                name=gene_set.name, genes=genes_in_overlap,
                warn_if_empty=False, representativeness=representativeness
            )
            for gene_set in self.gene_sets
            for genes_in_overlap in [gene_set.genes & genes]
            for representativeness in [len(genes_in_overlap) / len(gene_set.genes)]
            if representativeness >= min_representation
        })

    @property
    @lru_cache_mypy
    def all_genes(self) -> Set:
        all_genes = set()

        for gene_set in self.gene_sets:
            all_genes.update(gene_set.genes)

        return all_genes

    def to_frame(self, format='wide', include_metadata: bool = False) -> DataFrame:
        assert format in {'wide', 'long'}
        if format == 'wide':
            assert not include_metadata
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
        if format == 'long':
            return DataFrame(
                [
                    {
                        'name': gene_set.name,
                        'description': gene_set.description,
                        'genes': gene_set.genes,
                        **(
                            gene_set.metadata
                            if include_metadata else
                            {}
                        )
                    }
                    for gene_set in self.gene_sets
                ]
            ).set_index('name')

    @property
    @lru_cache_mypy
    def gene_sets_by_name(self) -> Dict[str, GeneSet]:
        by_name = {
            gene_set.name: gene_set
            for gene_set in self.gene_sets
        }
        assert len(self.gene_sets) == len(by_name)
        return by_name
    
    def to_rpy2(self):
        from rpy2.robjects.vectors import ListVector, StrVector

        data = {
            name: StrVector(list(gene_set.genes))
            for name, gene_set in self.gene_sets_by_name.items()
        }
        r_data = ListVector(data)
        return r_data
    
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
    def __init__(self, path, version='7.4'):
        self.path = Path(path)
        if not self.path.exists():
            raise ValueError(f'Could not find MSigDB: {self.path} does not exist')
        self.version = str(version)
        wildcard_path = str((self.path / f'*.v{self.version}.*.gmt').resolve())
        self.gene_sets = [
            self.parse_name(Path(p).name)
            for p in glob(wildcard_path)
        ]
        self.xml_path = None
        for candidate_dir in [self.path, self.path.parent]:
            candidate_path = candidate_dir / f'msigdb_v{version}.xml'
            if candidate_path.exists():
                self.xml_path = candidate_path
                break
        
    @staticmethod
    def add_metadata_from_xml(gene_sets: GeneSets, path: Union[Path, str]):
        tree = ElementTree.parse(str(path))
        root = tree.getroot()

        metadata_by_name = {
            child.attrib['STANDARD_NAME']: child.attrib
            for child in root.iter('GENESET')
        }
        for gene_set in gene_sets.gene_sets:
            gene_set.metadata = metadata_by_name[gene_set.name]

    def parse_name(self, name):
        parsed = re.match(rf'(?P<name>.*?)\.v{self.version}\.(?P<id_type>(entrez|symbols)).gmt', name)
        return parsed.groupdict()

    def resolve(self, gene_sets, id_type) -> Path:
        path = self.path / f'{gene_sets}.v{self.version}.{id_type}.gmt'
        if path.exists():
            return path
        else:
            raise ValueError(f'Unknown library: {path}!')

    def load(self, gene_sets: str, id_type: str) -> GeneSets:
        path = self.resolve(gene_sets=gene_sets, id_type=id_type)

        gene_sets = GeneSets.from_gmt(path, name=gene_sets)
    
        if self.xml_path:
            MolecularSignaturesDatabase.add_metadata_from_xml(gene_sets, self.xml_path)

        return gene_sets
