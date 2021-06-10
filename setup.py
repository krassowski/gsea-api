from setuptools import setup
from setuptools import find_packages


try:
    from pypandoc import convert

    def get_long_description(file_name):
        return convert(file_name, 'rst', 'md')

except ImportError:

    def get_long_description(file_name):
        with open(file_name) as f:
            return f.read()


if __name__ == '__main__':
    setup(
        name='gsea_api',
        packages=find_packages(),
        package_data={'gsea_api': ['py.typed']},
        version='0.3.4',
        license='MIT',
        description='Pandas API for Gene Set Enrichment Analysis in Python (GSEApy, cudaGSEA, GSEA)',
        long_description=get_long_description('README.md'),
        author='Michal Krassowski',
        author_email='krassowski.michal+pypi@gmail.com',
        url='https://github.com/krassowski/gsea-api',
        project_urls={
            'Bug Tracker': 'https://github.com/krassowski/gsea-api/issues',
            'Source Code': 'https://github.com/krassowski/gsea-api'
        },
        keywords=['gsea', 'gene', 'set', 'enrichment', 'cuda', 'pandas', 'api', 'GSEApy', 'cudaGSEA'],
        classifiers=[
            'Development Status :: 4 - Beta',
            'License :: OSI Approved :: MIT License',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX :: Linux',
            'Topic :: Utilities',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Software Development :: Libraries :: Python Modules',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7'
        ],
        install_requires=['numpy', 'pandas', 'jupyter_helpers'],
    )
