from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
    ext_modules=[Extension("cio", ["cio.c", "readcxsfile.c"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())

