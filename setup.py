from setuptools import setup

setup(
    name='hstools',
    version='2016.14.01',
    py_modules='hstools',
    install_requires=[
      'click',
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
      ],
    entry_points='''
      [console_scripts]
      hstools=run:cli
    ''',
    )
