



Examples pymie_core
------------------------------------


.. _example_pymie:

Example of Mie scattering sphere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    x   = 1.00                  # mie size parameter
    wl  = 8.15                  # wave length
    rad = wl*x/(2.0*np.pi)      # particle radius
    m   = complex(1.29,-0.0472) # refractive index

    l_angles = [0.0, 90./np.pi] # list of angles
    n_max = 0                   # if nmax = 0 the code estimates the cutoff of the series
                                # involved in the calculation.


    Qvec, albe, g_par, Cvec, Fmat, Smat = mie_parameters(wl, rad, m, l_angles, nmax)

        Qvec  = [Qsca, Qext, Qabs, Qbk, Qpr]
        Cvec  = [Csca, Cext, Cabs]
        Fmat  = [F11/Qsca, F12/Qsca, F33/Qsca, F34/Qsca]
        Smat  = [S1, S2]
        
    print('--------------------------------------')
    print(' m  (refraction index) = ',m   )
    print(' x  (mie parameter)    = ',x   )
    print(' wl (wave lenght)      =', wl  )
    print(' rad(radius)'          =', rad )
    print('--------------------------------------')

    print(' Qsca =', Qvec[0], '  Qext =', Qvec[1] )
    print(' Qabs =', Qvec[2], '  Qbk  =', Qvec[3] )
    print(' Qpr  =', Qvec[4])
    print(' Csca =', Cvec[0], '  Cext =', Cvec[1] )
    print(' Cabs =', Cvec[2], '  albe =', albe    )
    print(' gasy =', g_par )
    print('--------------------------------------')
    print('Angle values: ', langle)
    print('--------------------------------------')
    print(F11)
    print(F12)
    print(F33)
    print(F34)
    print('--------------------------------------')


.. _example_pymie_dsd:

Example of Mie scattering distribution spheres
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The calculation here shown is derived from the equation (2.48) in the reference 
:cite:`Hansen_1974`, see the pdf at [Hansen-1974](https://link.springer.com/article/10.1007/BF00168069).



.. code-block:: python

    def _mie_size_distribution(wl, langle, m, pdf, bins):
         """
         Example about how to estimate mie scattering properties of a distribution
         of spheres with sizes defined by a pdf on a set of bins.
         
         Integration method/function not full checked
         
         .. todo::
         
            Add an specific test of the implementation with scipy 
            
              - scipy.integrate import trapz
              - check also: 
                scipy.integrate.cumtrapz(y, x=None, dx=1.0, axis=-1, initial=None)
         
         """
         lQsca = []; lQext = []
         lQabs = []; lalbe = []
         nmax = 0
         for mix_rate, rad in zip(pdf, bins):
            out =  mie_parameters(wl, rad, m, langle, nmax)
            Qvec, albe, g_par, Cvec, Fmat, Smat = out
            Qsca, Qext, Qabs, Qbk, Qpr = tuple(Qvec)
            lQsca.append(Qsca*mix_rate)
            lQext.append(Qext*mix_rate)
            lQabs.append(Qabs*mix_rate)
            lalbe.append(albe*mix_rate)
            
         Qsca_dsd = int_trapz(lQsca, bins)
         Qext_dsd = int_trapz(lQext, bins)
         Qabs_dsd = int_trapz(lQabs, bins)
         albe_dsd = int_trapz(lalbe, bins)
         
         return  Qsca_dsd, Qext_dsd, Qabs_dsd, albe_dsd


