from setuptools import setup

setup(
    name='sarlacc',
    version='2015.11.10',
    py_modules='libsarlacc',
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
      'scikit-image',
      'scikit-learn',
      'scipy',
      ],
    entry_points='''
      [console_scripts]
      sarlacc=sarlacc:cli
    ''',
    )
