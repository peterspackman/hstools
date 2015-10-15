from setuptools import setup

setup(
    name='sarlacc',
    version='0.1',
    py_modules='libsarlacc',
    install_requires=[
      'Click',
      'numpy',
      'scipy',
      'matplotlib',
      'fastcluster',
      'h5py',
      'progressive',
      'periodictable',

      ],
    entry_points='''
      [console_scripts]
      sarlacc=sarlacc:cli
    ''',
    )
