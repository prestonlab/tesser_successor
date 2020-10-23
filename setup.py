import setuptools
from Cython.Build import cythonize

import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

def readme():
    with open('README.md') as f:
        return f.read()


setuptools.setup(
    name='tesser',
    version='0.1.0',
    description='Package for analyze Tesser community structure experiment.',
    long_description=readme(),
    long_description_content_type="text/markdown",
    license='GPLv3',
    url='http://github.com/prestonlab/tesser_successor',
    packages=setuptools.find_packages('tesser'),
    install_requires=[
        'numpy',
        'scipy',
        'networkx',
        'matplotlib',
        'pandas',
        'seaborn',
        'cython',
    ],
    ext_modules=cythonize(['tesser/csr.pyx',
                           'tesser/cfit.pyx']),
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.8',
    ]
)
