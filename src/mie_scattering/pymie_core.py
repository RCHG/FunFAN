# -*- coding: utf-8 -*-
"""
.. note::
   | **Program name**      : pymie_core.py
   | **Author**            : Ramiro Checa-Garcia
   | **email**             : r.checagarcia@gmail.com
   | **Purpose**           : Functions related with mie scattering calculations
   | **Revision History**  :
   |   - Nov-2014   v1.0     Contact: R.Checa
"""


__version__ = "PyMie v0.9 validated"

import sys

from scipy.special import jv
from scipy.special import yv
import numpy as np
from copy import deepcopy as deepc
import unittest


def _mie_angular(k, phi):
    """
    *Recurrent function to calculate Mie Angular functions*

    .. note::
       There are two functions :math:`\pi_{n}(\phi)` and :math:`\kappa_{n}(\phi)`
       The Mie angular functions are Legendre polynomials of :math:`cos(\phi)'.

    Args:
        k (integers): Order of the functions
        theta (float): Scattering Angle [degrees]

    Return:
        mie_p (float): Mie Function :math:`\pi_{k}` evaluated on :math:`\phi`.
        mie_t (float): Mie Function :math:`\kappa_{k}` evaluated on :math:`\phi`.
        
    .. hint::
    
       We have to keep in the equations l as float.

    .. note::
    
       This function is included for testing purposes because it is not efficient
       at all. For calculations of the Phase Matrix it is better to optimize it.
       
    """

    theta_rad = theta*np.pi/180.
    u = np.cos(theta_rad)
    n = float(k)
    l = n
    if k==0:
       mie_p = 0.0
       mie_t = 0.0
    elif k==1:
        mie_p = 1.0
        mie_t = u
    else:
        mie_p = ((2.*l-1)/(l-1))*u*mie_angular(n-1., theta)[0]- \
                (l/(l-1.))*mie_angular(n-2., theta)[0]
        mie_t = l*u*mie_p-(l+1)*mie_angular(n-1., theta)[0]

    return (mie_p, mie_t)


def _mie_angular_incre(k, phi, mie_1, mie_2):
    """
    *NON Recurrent function to calculate Mie Angular functions*

    .. note::
    
       This function was validated by checking
       the phase matrix global results on several cases.

    .. note::
    
       There are two functions :math:`\pi_{n}(\phi)` and :math:`\kappa_{n}(\phi)`
       The Mie angular functions are polynomial of :math:`cos(\phi)`, in fact,
       Legendre pol.

    Args:
        k (integers): Order of the functions
        theta (float): Scattering Angle [degrees]
        mie_1 (float): Mie Function :math:`\pi_{k-1}` evaluated on :math:`\phi`.
        mie_2 (float): Mie Function :math:`\kappa_{k-2}` evaluated on :math:`\phi`.
        
    Return:
        (tuple): mie-functions with index k for recurrency relations evaluated a :math:`\phi` 

          - mie_p (float): Mie Function :math:`\pi_{k}` evaluated on :math:`\phi`.
          - mie_t (float): Mie Function :math:`\kappa_{k}` evaluated on :math:`\phi`.
          - mie_1 (float): Mie Function :math:`\pi_{k-1}` evaluated on :math:`\phi`.
        
    .. hint::
    
       We have to keep in the equations l as float.

    .. note::
    
       The additional output is needed to avoid the recurrent function. Then it
       is possible to improve the performance.
       
    """

    theta_rad = phi*np.pi/180.
    u = np.cos(theta_rad)
    n = float(k)
    l = n
    if k==0:
       mie_p = 0.0
       mie_t = 0.0
    elif k==1:
        mie_p = 1.0
        mie_t = u
    else:
        mie_p = ((2.*l-1)/(l-1))*u*mie_1- (l/(l-1.))*mie_2
        mie_t = l*u*mie_p-(l+1)*mie_1

    return (mie_p, mie_t, mie_1)

def vec_mie_angular(nmax, phi):
    """
    *Function to calculate Mie Angular functions [1,NMAX]*
    
    Args:
        nmax (integers): Order of the functions
        phi (float): Scattering Angle [degrees]
        
    Return:
        (tuple): mie-functions with indices 0 to nmax evaluated a :math:`\phi` 

          - l_mie_p (list): Mie-function [0, nmax] :math:`\pi_{k}` evaluated on :math:`\phi`.
          - l_mie_t (list): Mie-function [0, nmax] :math:`\kappa_{k}` evaluated on :math:`\phi`.
    """
    l_mie_p = []
    l_mie_t = []
    mie_p1 = 0
    mie_p2 = 0

    for k in range(0, nmax+2):
        mie_p_k, mie_t_k, mie_p1 = _mie_angular_incre(k, phi, mie_p1, mie_p2)
        l_mie_p.append(mie_p_k)
        l_mie_t.append(mie_t_k)
        mie_p2 = mie_p1
        mie_p1 = mie_p_k

    return l_mie_p, l_mie_t


def mie_complex_coef(x, m, n, method='Bessel'):
    """
    *Function to calculate Mie Complex Coeffients an and bn*

    Args:
        x (float): mie size parameter
        m (complex): complex refraction index
        n (integer): order of an, bn
        
    Return:
        (tuple): mie complex coefficients :math:`(a_{n},b_{n})` of order n.

          - an (float): mie complex coefficient :math:`a_{n}`
          - bn (float): mie complex coefficient :math:`b_{n}`
        
    """
    jv_dn = jv(n-0.5, x)
    jv_up = jv(n+0.5, x)
    yv_dn = yv(n-0.5, x)
    yv_up = yv(n+0.5, x)
    arg_n = -n/(m*x)+jv(n-0.5, m*x)/jv(n+0.5, m*x)
    fac_na = arg_n/m+n/x
    fac_nb = m*arg_n+n/x

    an = (fac_na*jv_up-jv_dn)/complex(fac_na*jv_up-jv_dn, fac_na*yv_up-yv_dn)
    bn = (fac_nb*jv_up-jv_dn)/complex(fac_nb*jv_up-jv_dn, fac_nb*yv_up-yv_dn)

    return an, bn


def mie_parameters(wavelenght, rad, m, langle, nmax=0, eps='1e-8'):
    """
    *This function calculate all the physical parameters*

    .. todo::
       This function may be more optimized, the amounts an1, bn1 are also
       the an, bn of the next term. But here are always calculated an, bn,
       an1, bn1 for all terms. Also for each $\phi$ we have to calculate the
       Mie angular functions but it is not necessary calculate for all terms,
       technically only its necessary one time in the nval loop.

    Args:
        wavelenght (float): Wavelenght
        rad (float): Radius of the sphere
        m (complex): Complex Refraction Index
        langle (float): list of angles to calculate the phase matrix
        nmax (float): Max number of terms to be included on Mie calculations. If
                      nmax=0 the number of terms are calculated to provided an
                      estimated accuracy of eps.
        eps (float): By default it is 1e-8.
        
    Return:
        (tuple): mie scattering optical parameters. Note that Fmat and Smat elements are themselves lists.

          - Qvec  (list): [Qsca, Qext, Qabs, Qbk, Qpr]
          - albe  (float): single scattering albedo
          - g_par (float): asymmetry g parameter
          - Cvec  (list): [Csca, Cext, Cabs]
          - Fmat  (list): [F11/Qsca, F12/Qsca, F33/Qsca, F34/Qsca]
          - Smat  (list): [S1, S2]
        
    """

    wl = wavelenght
    x = 2.*np.pi*rad/wl
    if nmax == 0:
        nmax = int(x+4.3*x**(1./3.)+1)
    Qext = 0
    Qsca = 0
    Qbkp = 0
    newq = 0
    S1   = np.zeros(len(langle), dtype=complex)
    S2   = np.zeros(len(langle), dtype=complex)
    F11  = np.zeros(len(langle))
    F12  = np.zeros(len(langle))
    F33  = np.zeros(len(langle))
    F34  = np.zeros(len(langle))

    an=complex(0,0)
    bn=complex(0,0)

    amie_p = np.zeros((int(nmax)+3,len(langle)))
    amie_t = np.zeros((int(nmax)+3,len(langle)))

    for tindex, theta in enumerate(langle):
        amie_p[:,tindex], amie_t[:,tindex] = vec_mie_angular(int(nmax)+1, theta)

    for nval in range(1, int(nmax)+1):
        an, bn   = mie_complex_coef(x, m, nval)
        an1, bn1 = mie_complex_coef(x, m, nval+1) # Improve here performance.
        an1c = an1.conjugate()
        bn1c = bn1.conjugate()


        sum_ab = an + bn
        dif_ab = an - bn

        fac_n1 = float(nval)*(nval+2.)/(nval+1.)
        fac_n2 = (2.*float(nval)+1.)/(float(nval)*(float(nval)+1.))

        Qext = Qext + (2*float(nval)+1)*sum_ab.real
        Qsca = Qsca + (2*float(nval)+1)*(an.real**2+an.imag**2+bn.real**2+bn.imag**2)
        Qbkp = Qbkp + (2*float(nval)+1)*(-1)**nval*dif_ab
        newq = newq + fac_n1*((an*an1c+bn*bn1c).real) + fac_n2*((an*bn.conjugate()).real)

        for tindex, theta in enumerate(langle):
            
            S1[tindex]=S1[tindex]+(an*amie_p[nval,tindex]+bn*amie_t[nval,tindex])*fac_n2
            S2[tindex]=S2[tindex]+(an*amie_t[nval,tindex]+bn*amie_p[nval,tindex])*fac_n2

    for tind, theta in enumerate(langle):
        S1t = S1[tind]
        S2t = S2[tind]
        S1c = S1[tind].conjugate()
        S2c = S2[tind].conjugate()

        F11[tind]=(0.5*(+S1t.real**2 + S1t.imag**2 + S2t.real**2 + S2t.imag**2)).real
        F12[tind]=(0.5*(-S1t.real**2 - S1t.imag**2 + S2t.real**2 + S2t.imag**2)).real

        F33[tind]=(0.5*(S2c*S1t+S2t*S1c)).real
        F34[tind]=(0.5*(S1t*S2c-S2t*S1c)*complex(0,1)).real


    fac   = 2./x**2
    Qbk   = 0.5*fac*(Qbkp.real**2+Qbkp.imag**2)
    Qsca  = fac*Qsca
    Qext  = fac*Qext
    gQsca = 2.*fac*newq
    g_par = gQsca/Qsca
    albe  = Qsca/Qext
    Qabs  = Qext-Qsca
    Qpr   = Qext-gQsca
    Csca  = Qsca*np.pi*rad**2
    Cext  = Qext*np.pi*rad**2
    Cabs  = Qabs*np.pi*rad**2

    Qvec  = [Qsca, Qext, Qabs, Qbk, Qpr]
    Cvec  = [Csca, Cext, Cabs]
    Fmat  = [F11/Qsca, F12/Qsca, F33/Qsca, F34/Qsca]
    Smat  = [S1, S2]

    return Qvec, albe, g_par, Cvec, Fmat, Smat


def _mie_size_distribution(wl, langle, m, pdf, bins):
     """
     Example about how to estimate mie scattering properties of a distribution
     of spheres with sizes defined by a pdf on a set of bins.

     Note about integration:
     
     scipy.integrate.cumtrapz(y, x=None, dx=1.0, axis=-1, initial=None)

     """
     from scipy.integrate import trapz as int_trapz

     lQsca = []; lQext = []
     lQabs = []; lalbe = []
     nmax = 0
     for mix_rate, rad in zip(pdf, bins):
        Qvec, albe, g_par, Cvec, Fmat, Smat =  mie_parameters(wl, rad, m, langle, nmax)
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



class Test(unittest.TestCase):

    def test_mie_params_noabs(self):
        x  = 4.9646
        wl = 0.6328
        rd = wl*x/(2.0*np.pi)
        m  = complex(1.5,0)
        langle= [-0.6*180.0/np.pi]
        Qvec, albe, g_par, Cvec, Fmat, Smat  = mie_parameters(wl, rd, m, langle, nmax=0, eps='1e-8')
        # Qvec  = [Qsca, Qext, Qabs, Qbk, Qpr]
        self.assertAlmostEqual( Qvec[0]  ,  3.89618111 , 5 )
        self.assertAlmostEqual( Qvec[1]  ,  3.89618111 , 5 )
        self.assertAlmostEqual( Qvec[2]  ,  0.0        , 5 )
        self.assertAlmostEqual( Qvec[3]  ,  1.94290317 , 5 )
        self.assertAlmostEqual( Qvec[4]  ,  1.13903412 , 5 )
        self.assertAlmostEqual( g_par    ,  0.70765370 , 5 )

    def test_mie_params_abs(self):
        x  = 4.9646
        wl = 0.6328
        rd = wl*x/(2.0*np.pi)
        m  = complex(1.5,-0.01)
        langle= [0.6*180.0/np.pi]
        Qvec, albe, g_par, Cvec, Fmat, Smat  = mie_parameters(wl, rd, m, langle, nmax=0, eps='1e-8')
        # Qvec  = [Qsca, Qext, Qabs, Qbk, Qpr]
        self.assertAlmostEqual( Qvec[0]  ,  4.32425720 , 5 )
        self.assertAlmostEqual( Qvec[1]  ,  3.99370611 , 5 )
        self.assertAlmostEqual( Qvec[2]  , -0.33055108 , 5 )
        self.assertAlmostEqual( Qvec[3]  ,  2.83593856 , 5 )
        self.assertAlmostEqual( Qvec[4]  ,  1.05595487 , 5 )
        self.assertAlmostEqual( g_par    ,  0.67936551 , 5 )
        
if __name__ == '__main__':
    unittest.main()




