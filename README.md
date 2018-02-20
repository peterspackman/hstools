# hstools
This repository contains a set of tools for spherical harmonic transforms and the
description of isosurface/surface properties. This library is available under the
conditions of the GPLv3 (see LICENSE).

Installation should simply be a matter of:
```
pip install git+https://github.com/peterspackman/hstools
```

*Note that this is experimental software, and as such may require more effort/understanding 
to use. While I aim to make the library as user-friendly as I can, it's a work in progress.*


## Additional information
This software can optionally make use of the [shtns](https://users.isterre.fr/nschaeff/SHTns/) 
library, which I recommend using. If you do so, please additionally cite this paper:
```
@article {shtns,
  author = {Schaeffer, Nathanael},
  title = {Efficient spherical harmonic transforms aimed at
  pseudospectral numerical simulations},
  journal = {Geochemistry, Geophysics, Geosystems},
  doi = {10.1002/ggge.20071},
  volume = {14}, number = {3}, pages = {751--758},
  year = {2013},
}
```
