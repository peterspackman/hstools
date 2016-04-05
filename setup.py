from setuptools import setup

setup(
    name='hstools',
    description=("Module and CLI for processing HDF5 HS files"),
    version='2016.04.04-2',
    packages=find_packages(),
    install_requires=[
      'fastcluster',
      'h5py',
      'hdbscan',
      'matplotlib',
      'numpy',
      'periodictable',
      'plyfile',
      'progressive',
      'scikit-learn',
      'scipy',
      'tqdm',
      ],
    entry_points= {
      'console_scripts': ['hstools=hstools.command_line:main'],
    }
)
