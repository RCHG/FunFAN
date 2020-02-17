

.. _xarray: http://xarray.pydata.org
.. _pandas: http://pandas.pydata.org
.. _netCDF: http://www.unidata.ucar.edu/software/netcdf
.. _matplotlib: https://matplotlib.org/
.. _cartopy: https://scitools.org.uk/cartopy/docs/latest/



Set of Functions
================

Particle size distributions (aersols_psd) 
------------------------------------------

.. note::

    - *Author* R. Checa-Garcia
    - *License* GPLv3 (GNU General Public License version 3)
    - *Date* 2020-01-01 (last update see github)
    - *Other requirements* xarray_ , matplotlib_


Set of functions to get **statistics** information of the modal particle size distributions,
create **plots** (using matplotlib_) with the modal distributions, and perform a **binning** of typical aerosols diagnostics.
To manage the netCDF files for binning methods, it is used the library xarray_

 - Log-normal distribution
 - Binning of (time, lat, lon) fields. Tests for tendency fluxes (emissions, dry and wet deposition)



.. seealso::
   
   - Examples :ref:`ln_example`
   - Examples :ref:`bins_emi_example`
   - Tests
   
         - If you run ``python3 lognormal.py`` the code perform several unittests.
         - The code for bins has been tests, the example is at the code.
         - There are few private methods (not in documentation but in the code) with examples/tests.

