"""
.. note::
   | **Program name**      : lognormal.py
   | **Author**            : Ramiro Checa-Garcia
   | **email**             : rcheca@lsce.ipsl.fr
   | **Purpose**           : Functions related with lognormal aerosol distributions.
   | **Revision History**  :
   |   - Jan-2020   v1.0     Initial version. Contact: R.Checa
"""

import xarray as xr
import numpy as np
from scipy.special import erf
import unittest

def Xpdf_logD(D, lnmu, lnsigma):
    '''
    Function for the pdf (probability density function) 
     :math:`x^{*}(D)` such as :math:`dX(D)=x^{*}(D; ln \mu, ln \sigma)dlnD`

    Args:
        D (float): diameter value
        lnmu (float): mean value as :math:`ln(\mu)`
        lnsigma (float): standard deviation as :math:`ln(\sigma)`

    Returns:
        pdf (float):Pprobability at diameter D (but note about dlnD, so we have actually a gaussian distribution).
        
    Comment:
        This function is actually very general as:

         - :math:`dN(D)=n^{*}(D; ln \mu_{nmd}, ln \sigma)dlnD`
         - :math:`dA(D)=a^{*}(D; ln \mu_{amd}, ln \sigma)dlnD`
         - :math:`dV(D)=v^{*}(D; ln \mu_{mmd}, ln \sigma)dlnD`

        where

         - :math:`\mu_{nmd}` (number median diameter)
         - :math:`\mu_{amd}` (area   median diameter)
         - :math:`\mu_{mmd}` (volume median diameter)
        
    '''

    pdf = np.exp(-(np.log(D) - lnmu)**2 / (2.0 * lnsigma**2)) / (lnsigma * np.sqrt(2.0 * np.pi))
    
    return pdf

def Npdf_D(D, lnmu, lnsigma):
    '''
    Function for the pdf (probability density function) 
     :math:`x^{*}(D)` such as :math:`dX(D)=x^{*}(D; ln \mu, ln \sigma)dD`

    Args:
        D (float): diameter value
        lnmu (float): mean value as :math:`ln(\mu)`
        lnsigma (float): standard deviation as :math:`ln(\sigma)`

    Returns:
        pdf (float): Probability at diameter D (note about dD)
        
    '''
    
    pdf = D**(-1)*np.exp(-(np.log(D) - lnmu)**2 / (2.0 * lnsigma**2)) / (lnsigma * np.sqrt(2 * np.pi))
    return pdf


def Xcdf_D(D, lnmu, lnsigma):
    '''
    Function for the cdf (cumulative distribution function) 
    :math:`X^{*}(D)` such as :math:`X(D)=\int_{0}^{D}x^{*}(y; ln \mu, ln \sigma)dy`

    Args:
        D (float): diameter value
        lnmu (float): mean value as :math:`ln(\mu)`
        lnsigma (float): standard deviation as :math:`ln(\sigma)`

    Returns:
        cdf (float): Cumulative Probability Function from 0 to D.
        
    '''
    
    cdf = 0.5*(1+erf((np.log(D) - lnmu) / (np.sqrt(2.0)*lnsigma)))

    return cdf


def plot_multimodal(ax0, a_sigma, a_mmd, a_names,  a_fact, nmodes,
                    dmax=150.0, dmin=0.001, dinc=0.0005,
                    rho=2650, nconc=7e6, ndist='v'):
    """
    This function add a plot of a multimodal distribution ot the axes ax0,
    the modes added are those defined by a_sigma and a_mmd but only using
    the first nmodes in the array.
    
    Args:
        ax0 (matplotlib axes): A matplotlib axes where draw the distribution
        a_sigma (array): Array with the values of sigma of a multimodal lognormal.
        a_mmd (array): Array with the values of the mass median diameter of all modes.
        a_names (array): Array or list with the names of the modes.
        a_fact (array): Array with fractions for each mode.
        nmodes (array): Index (in the arrays above) of modes to use.
        dmax (float, optional): Max diameter on plot. Defaults to 150.0
        dmin (float, optional): Min diameter on plot. Defaults to 0.001
        dinc (float, optional): resolution in diameter. Defaults to 0.0005
        rho (float, optional): Density of the particles. Defaults to 2650
        nconc (float, optional): Number concentration. Defaults to 7e6
        ndist (string, optical): Name of distribution. Defaults to v.

    Returns:
        ax0 (matplotlib axes): return the axes introduced and changed by the function.

    """
    
    a_lnsigma = np.log(a_sigma)
    x = np.arange(dmin, dmax, dinc)
    pltsum_V = np.zeros_like(x)

    for ii in nmodes:
        ynew_V = a_fact[ii]*Xpdf_logD(x, np.log(a_mmd[ii]), a_lnsigma[ii])
        ax0.plot(x, ynew_V, lw=1.5, ls='--', label=a_names[ii])  # each mode
        pltsum_V = pltsum_V + ynew_V

    ax0.plot(x, pltsum_V, lw=3.0) # multimodal
    ax0.set_ylabel(r'$\frac{d'+ndist+'(D)}{dlnD}$', fontsize=14, rotation='horizontal', labelpad=14)

    return ax0


    
    
def stats_lognormal_mode(sigma, mmd, name, rhop=2650.0, n_conc=7.e6, show=True):
    '''
    Estimates, and optionally shows, the general info of an aerosol lognormal mode.

    Args:
        sigma (float): Standard deviation of mode.
        mmd (float): Mass median diameter.
        rhop (float, optional): Density particles in mode. Defaults to 2650.0.
        n_conc (float, optional): Number concentration. Defaults to 7e6
        show (bool, optional): True to print in stdout. Defaults to True.

    Returns:
        out (list): list with several statistics asociated to lognormal distr.

    .. testcode::

       # =======================================================================
       # Typical output 
       #
       # > stats_lognormal_mode(1.8, 1e-6, 'D1', rhop=2650.0, n_conc=7.e6)

       =========================================================================
       *sigma  = 1.8000e+00  | *mmd   = 1.0000e-06 m   |=> nmd    = 3.5470e-07 m 
                             |=> smd  = 7.0787e-07 m   |   smd!   = 7.0787e-07 m 
        reff_n = 8.4135e-07  | reff_v = 5.0000e-07 m   |   reff_m = 5.0000e-07 m 
        mmeand = 1.1886e-06  | nmoded = 2.5108e-07 m   |   d_avgm = 5.9557e-07 m 
       *n_conc = 7.0000e+06  | a_conc = 4.6456e-06 m   |   v_conc = 7.7427e-13 m 
        m_conc = 2.0518e-09   
    
       *before word means pre-defined, ! after word means 2nd method 
       all others without symbols mean calculated
       =========================================================================
    '''

    volume  = (1.0/6.0)*np.pi*((mmd)**3)
    surface = 4.0*np.pi*((0.5*mmd)**2)
    reff_va = 3.0*(1.0/6.0)*np.pi*((mmd)**3)/(4.0*np.pi*((0.5*mmd)**2))
    lnsigma = np.log(sigma)

    nmd = mmd*np.exp(-3.0*lnsigma**2)
    smd = mmd*np.exp(-lnsigma**2)
    smd_bis = nmd*np.exp(2.0*lnsigma**2)  # same than smd but other way

    reff_n= nmd*np.exp(2.5*lnsigma**2)    # eff radius of n(D)
    reff_v= reff_va                       # eff radius of v(D)
    reff_m= mmd/2                         # eff radius of m(D)
    
    mmeand = nmd*np.exp( 3.5*lnsigma**2)  # mass-weighted mean diameter   -> right this is the mass-weighted mean diameter
                                          #                                  (different of mass median diameter = mmd)
                                          # formally the mass mean diameter is what is named below the diameter of average mass
    nmoded = nmd*np.exp(-1.0*lnsigma**2)  # number mode diameter
    d_avgm = nmd*np.exp( 1.5*lnsigma**2)  # diameter of average mass

    #n_conc = (surface*(nmd**(-2))/np.pi)*np.exp(1.5*lnsigma**2)
    a_conc = (nmd**2)*n_conc*np.pi*np.exp(1.5*lnsigma**2)
    m_conc = (nmd**3)*(n_conc/6.0*np.pi)*rhop*np.exp(4.5*lnsigma**2)
    v_conc = (nmd**3)*(n_conc/6.0*np.pi)*np.exp(4.5*lnsigma**2) 

    
    if show==True:
        print('\n =========================================================================')
        print(' *sigma  = %4.4e  | *mmd   = %4.4e m   |=> nmd    = %4.4e m ' % (sigma, mmd, nmd) )
        print('                       |=> smd  = %4.4e m   |   smd!   = %4.4e m ' % (smd, smd_bis) )
        print('  reff_n = %4.4e  | reff_v = %4.4e m   |   reff_m = %4.4e m ' % (reff_n, reff_v, reff_m) )
        print('  mmeand = %4.4e  | nmoded = %4.4e m   |   d_avgm = %4.4e m ' % (mmeand, nmoded, d_avgm) )
        print(' *n_conc = %4.4e  | a_conc = %4.4e m   |   v_conc = %4.4e m ' % (n_conc, a_conc, v_conc) )
        print('  m_conc = %4.4e   ' % (m_conc) )
        print(' \n ... info ...')
        print(' * before word means pre-defined, ! after word means 2nd method ')
        print(' all others without symbols mean calculated')
        print(' =========================================================================\n')

    out = [sigma, mmd, nmd, smd, smd_bis, reff_n, reff_v, reff_m,
           mmeand, nmoded, d_avgm, a_conc, v_conc, m_conc]
    
    return out


def bin_fractions_lognormal(xmd, sigma, bins):
    """
    Function to estimate the bins fractions for given :math:`\mu_{md}` and :math:`\sigma`.
    Note that :math:`\mu_{md}` is usually :math:`\mu_{nmd}`, :math:`\mu_{amd}` or :math:`\mu_{mmd}`    
    It assumed a unit consistency between xmd and the bins values. 

    Args:
        xmd (float): x median diameter
        sigma (float): :math:`\sigma` (not :math:`ln \sigma`).
        bins (list): List of bin limits (no central values).

    Returns:
        lfractions (list): list with bins fraction between bin limits.

    """
    lnsigma = np.log(sigma)
    lfractions = []
    for x0,x1 in zip(bins[0:-1],bins[1::]):
        fraction = Xcdf_D(x1, np.log(xmd), lnsigma)-Xcdf_D(x0,  np.log(xmd), lnsigma)
        lfractions.append(fraction)
    
    return lfractions

class lnmode(object):
    """ 
    Object with key information of an atmospheric particle lognormal mode
    It is based on a dictionary used to build de object, this dictionary
    has some fields as mandatory:

        - ``sigma``  :math:`\sigma` of lognormal mode (not :math:`ln \sigma`!!)
        - ``n_conc`` number concentration
        - ``name`` name of mode
        - ``rhop`` density particle

    and one of these two:

        - ``mmd``  mass median diameter
        - ``nmd``  number median diameter
            
    Args:
        ini_dict (dict): Dictionary with the main information to be used in the lognormal mode. It is mandatory that this dataset has the information indicated above.


    .. testcode::

        # ======= EXAMPLE of how to use (piece of code, not in REPL)
        #
        # (1) Create/Define the object here a lnmode named lnD1
        lnD1 = lnmode({'sigma': 1.8, 'mmd':1.e-6, 'n_conc':7e-6, 'rhop':2650.0, 'name':'D1'})

        # (2) How the statistics information of this mode
        lnD1.show()

        # (3) Estimate the bin fractions for that distribution for mass distribution 
        print(lnD1.bin_fractions([0.001e-6,1.0e-6,2.e-6,3.e-6]))

        # (4) Estimate the bin fractions for that distribution for number distribution
        print(lnD1.bin_fractions([0.001e-6,1.0e-6,2.e-6,3.e-6]), kind='nmd')

    """

    def __init__(self, ini_dict):
        self.sigma = ini_dict['sigma']
        self.n_conc = ini_dict['n_conc']
        self.name = ini_dict['name']
        self.rhop = ini_dict['rhop']
        self.lnsigma = np.log(self.sigma)
        if 'mmd' in ini_dict.keys():
            self.mmd = ini_dict['mmd']
            self.nmd = self.mmd*np.exp(-3.0*self.lnsigma**2)
            self.smd = self.mmd*np.exp(-self.lnsigma**2)

        elif 'nmd' in ini_dict.keys():
            self.nmd = ini_dict['nmd']
            self.mmd = self.nmd*np.exp(3.0*self.lnsigma**2)
            self.smd = self.mmd*np.exp(-self.lnsigma**2)
        else:
            print(' Problem here')

    def show(self):
        """
        Shows the main properties of the lognormal mode
        """
        stats_lognormal_mode(self.sigma, self.mmd, self.name, self.rhop, self.n_conc, show=True)

    #def plot(self, plotname):


    def bin_fractions(self, bins, kind='mmd'):
        """
        Calculate the fraction of the distribution of each bin. We have an
        object with mmd, nmd etc. So we have to be explicit about for which
        distribution we want to calculate bin fractions (mmd, nmd etc).
        
        NOTE: we have to be consistent between bins units and mmd/nmd units.

        Args:
            bins (list/array): List with all bin limits (no bins center values)
            kind (str, optional): mmd or nmd. Defaults to 'mmd'.

        Returns:
            bin_frac (list): Freq/fraction per bin

        """

        if kind=='mmd':
            lnd = np.log(self.mmd)
        elif kind=='smd':
            lnd = np.log(self.smd)
        elif kind=='nmd':
            lnd = np.log(self.nmd)

        bin_frac = [ Xcdf_D(x1, lnd, self.lnsigma)- Xcdf_D(x0, lnd, self.lnsigma)
                    for x0,x1 in zip(bins[0:-1],bins[1::])]

        return bin_frac


def _check_pdf_cdf(dmin, dmax, lnmu, lnsigma, str1='n', show=True):
    '''
    This function has been created for testing purposes. It was designed to 
    compare values of an excel sheet with lookup tables of modal distributions.
    We used to ascertain that our functions pdf, cdf are consistent.

    Args:
        dmin (float): min diameter in interval.
        dmax (float): max diameter in interval.
        lnmu (float): mean value as mu
        lnsigma (float): standard deviation as sigma
        str1 (str, optional): n for number distr. v for volume distr. Defaults to 'n'.
        show (bool, optional): True to print in stdout. Defaults to True.

    Returns:
        None.

    #print(' diameter ----- dV/dlnD ----- dV/dD ')
    # COMPARISON WITH EXCEL FILES:
    #
    # Npdf_D(3.0, 2.5, 1.02) is consistent with LOGNORM.DIST(3.0, 2.5, 1.02, FALSE)
    #
    # Therefore:
    #
    # check_pdf_cdf(dmin, dmax, np.log(nmd*1e6), lnsigma, str1='n')
    # check_pdf_cdf(dmin, dmax, np.log(mmd*1e6), lnsigma, str1='v')
    #
    # reproduce the results of the excel of Michal Schulz but with python.
    # In the excel the PDF is estimated from the CDF, with the custom_pdf style, but there
    # we show that the this custom_pdf for a small interval is a value between pdf(dmin) and
    # pdf(dmax), (remember here the mean value theorem that there is a d_inside with the true value.
    # so the custom_pdf is a non accurate version of the actual pdf.
    #

    for pair in [(0.020, 0.021), (0.206, 0.220), (4.244, 4.539)]:
        dmin = pair[0]
        dmax = pair[1]
        an, bn = check_pdf_cdf(dmin, dmax, np.log(nmd*1e6), lnsigma, str1='n', show=show)
        av, bv = check_pdf_cdf(dmin, dmax, np.log(mmd*1e6), lnsigma, str1='v', show=show)

    #######
    (0.020, 0.021) -> an = 4.85e-6 || av = 1.93e-10 || bn = 2.40e-4 (values approx. based on cdf)
                           5.31e-6         2.17e-10         2.59e-4 (values using pdf)
    (0.206, 0.220) -> an = 4.66e-1 || av = 2.14e-2  || bn = 2.19
                           4.66e-1         2.13e-2          2.19    (values using pdf)
    (4.244, 4.539) -> an = 7.22e-5 || av = 2.87e-2  || bn = 1.65e-5
                           7.22e-5         2.87e-3          1.64e-5 (values using pdf)

    For smaller values than 1.e-6 there are more differences due to approx.
    
    '''
    
    anlog = Xcdf_D(dmin, lnmu, lnsigma)
    bnlog = Xcdf_D(dmax, lnmu, lnsigma)
    cnlog = (bnlog-anlog)/(np.log(dmax)-np.log(dmin))

    anD = Xcdf_D(dmin, lnmu, lnsigma)
    bnD = Xcdf_D(dmax, lnmu, lnsigma)
    cnD = (bnD-anD)/(dmax-dmin)

    apdf= Npdf_D(dmin, lnmu, lnsigma)
    bpdf= Npdf_D(dmax, lnmu, lnsigma)

    if show==True:
        if str1=='n':
            print(      ' Dmin = %4.3f  Dmax = %4.3f     dn/dlnD = %4.4e      dn/dD  = %4.4e ' % (dmin, dmax, cnlog, cnD))
            print(      ' Dmin = %4.3f  Dmax = %4.3f     n(dmin) = %4.4e      n(dmax)= %4.4e    dn = %4.4e' % (dmin, dmax, apdf, bpdf, 0.5*(bpdf+apdf)*(dmax-dmin)))
    
        elif str1=='v':
            print(      ' Dmin = %4.3f  Dmax = %4.3f     dv/dlnD = %4.4e      dv/dD  = %4.4e ' % (dmin, dmax, cnlog, cnD))
            print(      ' Dmin = %4.3f  Dmax = %4.3f     v(dmin) = %4.4e      v(dmax)= %4.4e    dv = %4.4e' % (dmin, dmax, apdf, bpdf, 0.5*(bpdf+apdf)*(dmax-dmin)))
        else:
            print('--')
        print('\n')

    return cnlog, cnD

class Test(unittest.TestCase):

    def _test_info_lognormal(self):
        test = stats_lognormal_mode(1.8, 1e-6, 'D1', rhop=2650.0, n_conc=7.e6)
        (sigma, mmd, nmd, smd, smd_bis, reff_n, reff_v, reff_m, mmeand, nmoded, d_avgm, a_conc, v_conc, m_conc) = test
        self.assertEqual(       sigma  ,  1.8 )
        self.assertEqual(       mmd    ,  1e-6 )
        self.assertAlmostEqual( nmd    ,  3.55e-7 ,9 )
        self.assertAlmostEqual( smd    ,  7.08e-7 ,9)
        self.assertAlmostEqual( smd_bis,  7.08e-7 ,9)
        self.assertAlmostEqual( reff_n ,  8.41e-7 ,9)
        self.assertAlmostEqual( reff_v ,  5.0e-7  ,9)
        self.assertAlmostEqual( reff_m ,  5.0e-7  ,9)
        self.assertAlmostEqual( mmeand ,  1.19e-6 ,8)
        self.assertAlmostEqual( nmoded ,  2.51e-7 ,9)
        self.assertAlmostEqual( d_avgm ,  5.96e-7 ,9)
        self.assertAlmostEqual( a_conc ,  4.65e-6 ,8)
        self.assertAlmostEqual( m_conc ,  2.05e-9 ,11)
        self.assertAlmostEqual( v_conc ,  7.74e-13,14)

if __name__ == '__main__':
    unittest.main()
