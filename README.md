## Hirshfeld surface binary

```
git submodule update --init --recursive
mkdir build && cd build
cmake ..
make -j4
./hs -i input.cif -o output.h5 -b basis_sets -l lmax
```
