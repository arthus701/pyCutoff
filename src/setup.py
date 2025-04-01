import os
import codecs
import setuptools


# https://packaging.python.org/guides/single-sourcing-package-version/
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


long_description = ''

name = 'pyCutoff'
version = get_version('pyCutoff/__init__.py')
description = '''Cutoff rigidities in python'''
copyright = '2025 Helmholtz Centre Potsdam GFZ, ' \
    + 'German Research Centre for Geosciences, Potsdam, Germany'


setuptools.setup(
    name=name,
    version=version,
    author='Schanner, M. A.',
    author_email='arthus@gfz.de',
    packages=['pyCutoff'],
    license='GPL v3',
    description=description,
    long_description=long_description,
    install_requires=[
        'numpy>=1.18',
        'tqdm',
        'paleokalmag',
        'pymagglobal',
        'pyproj',
        ],
)
