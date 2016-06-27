## Hirshfeld surface binary

```
git submodule update --init --recursive
mkdir build && cd build
cmake ..
make -j4
./hs -i input.cif -o output.h5 -b basis_sets -l lmax
```
### HSTools

The aim of this project is to provide a simple interface
for processing data from Hirshfeld surface calculations
and their associated spherical harmonic shape descriptors.

This project is written for Python 3.

# Installation
Installing should be as follows:
```
pip3 install git+git://github.com/peterspackman/hstools
```

If you're using conda you may need to install the dependencies
before trying to install this software

# Dependencies
* fastcluster
* hdbscan
* h5py

* scikit-learn
* scipy
* matplotlib
* numpy

* periodictable
* plyfile
* tqdm
