from setuptools import setup, find_packages

setup(name='hstools',
      version='0.1.1',
      description='Tools for analysing hirshfeld surfaces',
      url='http://github.com/peterspackman/hstools',
      author='Peter Spackman',
      author_email = 'peterspackman@fastmail.com',
      keywords = ['crystallography', 'hirshfeld surface', 'spherical harmonic transform'],
      classifiers = [
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Development Status :: 2 - Pre-Alpha',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Software Development :: Libraries :: Python Modules',
          ],
      license='GPLv3',
      packages=find_packages(),
      package_data={
         'hstools': ['*.h5']  # include all templates
      },
      install_requires=['numpy', 'h5py', 'sbf', 'scipy', 'trimesh'],
      entry_points={
          'console_scripts': [
              'hsdecompose = hstools.decompose:main',
          ]
      },
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose']
)
