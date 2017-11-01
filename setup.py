from setuptools import setup, find_packages

setup(name='hstools',
      version='0.1.7',
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
         'hstools': ['*.sbf'],  # include all sbf files,
         'hstools': ['*.bin']  # and all numpy array files
      },
      install_requires=['numpy', 'sbf', 'scipy', 'pandas'],
      entry_points={
          'console_scripts': [
              'hstools-describe = hstools.decompose:main',
              'hstools-search = hstools.search:main',
              'hstools-fakecif = hstools.fakecif:main',
              'hstools-visualize= hstools.visualize:main',
          ]
      },
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose']
)
