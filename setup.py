from setuptools import setup, find_packages

setup(
    name='hstools',
    version='2016.04.05',
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
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
    package_data = {
        'data': ['*.h5'],
    },
    author='Peter Spackman',
    author_email='20265845@student.uwa.edu.au',
    license='BSD',
    description=("Module and CLI for processing HDF5 HS files"),
    entry_points= {
      'console_scripts': ['hstools=hstools.command_line:main'],
    }
)
