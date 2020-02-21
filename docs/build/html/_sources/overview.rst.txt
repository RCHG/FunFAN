

Overview
====================

This is set of modules and **fun**\ ctions in **P**\ ython used for the
analysis of **A**\ erosols in climate models. Rather than a full package 
that constrains to the user to adopt an specific work-flow, the
idea is to provide functions that can be implemented in other projects
with minimal changes. The common shared set of functions would help
to improve the reproducibility of our scientific results.

The package or set of functions is organized in folders where the
specific goal is described, as well as its requirements in terms of
external libraries.


General requirements
--------------------

Common set requirements to all different functions:

* ``Python 3``  it is based on python 3 (not python 2 syntax).
* ``numpy`` for manage of numerical and masked arrays.
* ``scipy`` for special functions and statistics methods.
* ``os``    library to manage file-system access with python.


Acknowledgements
----------------

- *Jasper Kok* for fostering discussions about binning.
- *Yves Balkanski* for providing excel files to test several calculations (of lognormals)
- *Andre Butz* for discussion about Mie scattering in the context of remote sensing of GHGs.



.. _xarray: http://xarray.pydata.org
.. _pandas: http://pandas.pydata.org
.. _netCDF: http://www.unidata.ucar.edu/software/netcdf
.. _matplotlib: https://matplotlib.org/
.. _cartopy: https://scitools.org.uk/cartopy/docs/latest/

