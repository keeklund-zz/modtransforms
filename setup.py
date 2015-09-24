""" """
from setuptools import setup, find_packages

__version__ = '0.0.3'

# setup packages to be installed
packages = find_packages(exclude=['*.tests',])

package_data = {}

# include useful and/or necessary scripts
scripts = ['bin/modtransform_simple',
        'bin/modtransform_seq',]

# python packages required for install
install_requires = [
    "argparse >= 1.3.0",
    "pysam >= 0.8.2.1",
]

setup_args = {
    "name"             : "modtransforms",
    "version"          : __version__,
    "author"           : "Karl Eklund",
    "author_email"     : "keklund@ad.unc.edu",
    "description"      : "Update various file formats using MOD File",
    "packages"         : packages,
    "package_data"     : package_data,
    "scripts"          : scripts,
    "license"          : open("LICENSE").read(),
    "url"              : "",
    "install_requires" : install_requires,
}

setup(**setup_args)
