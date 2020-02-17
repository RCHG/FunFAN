"""
.. note::
   | **Program name**      : binning.py
   | **Author**            : Ramiro Checa-Garcia
   | **email**             : rcheca@lsce.ipsl.fr
   | **Purpose**           : Functions to perform a binning of a modal based diagnostics
   | **Revision History**  :
   |   - Jan-2020   v1.0     Initial version. Contact: R.Checa
"""

import xarray as xr
import numpy as np
from scipy.special import erf
import os.path

import aerosols_psd.lognormal as ln 


def total_tendency(tendency_data, area):
    """
    Estimate the total tendency for a year
    assuming a year with 365 days. Inputs are in SI, and outputs are
    Tg per year.


    Args:
        tendency_data (datarray): Datarray with a tendency field with monthly resolution, and SI units.
        area (datarray):          Datarray with the area per grid cell.

    Returns:
        tendency (float):         Total tendency per year in Tg.

    """
    
    tendency_field = tendency_data*area
    tendency_month = tendency_field.sum(['lat','lon'])
    daysmon        = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    sec_month      = np.array([k * 86400. for k in daysmon])
    tendency_total = np.array([tendency_month[imon]*sec_month[imon] for imon in range(12)])
    
    return np.round(tendency_total.sum()/1.e9,8)

def _bins_by_model():
    """
    Dictionaries with typical bins of several models.
    Note that some of them has been extended with one or more additional bins.

    Returns:
        dic_bins (dic): prescribed bins of several models.

    """
    dic_bins = {
        # Original
        'CNRM-dust':   [0.0001, 0.1,    0.2, 0.5,   1.0, 2.5, 10.0, 100.0], 
        'UKESM-dust':  [0.0001, 0.0722, 0.2, 0.722, 2.0, 7.22, 20.0, 72.2],
        'GOCART-dust': [0.0001, 0.2, 2.0, 3.6, 6.0, 12.0, 20.0],
         
        # Extended
        'CNRM-dust-ext':   [0.0001, 0.1,    0.2, 0.5,   1.0, 2.5, 10.0, 100.0, 200.0], 
        'UKESM-dust-ext':  [0.0001, 0.0722, 0.2, 0.722, 2.0, 7.22, 20.0, 72.2, 200.0],
        'GOCART-dust-ext': [0.0001, 0.2, 2.0, 3.6, 6.0, 12.0, 20.0, 30.0, 40.0,
                            100.0, 200.0, 300.0]
               }

    return dic_bins


def create_netcdf_2D_bins(dic_fractions, dic_varname, dic_ncnames, modes, new_bins,
                           dic_files,  newvarname, ftest=None, save=True,
                           test_info=None):
    """
    Creates a netcdf of a given (lat,lon) diagnostic with values over a set of bins, given the
    netcdfs over a set of modes for that diagnostics. The inputs are therefore a set
    of netcdfs, one per mode, for a variable X(time, lat, lon) the output is a single netcdf with 
    X(time, lat, lon, bins).

    For the calculation we need the expected fractions per bin for each mode. That can be calculated
    with other functions provided in aerosols_psd set of functions.

    With test_info the function can perform test of consistency.

    Args:
        dic_fractions (dict):  Dictionary for fractions per mode for new_bins.
        dic_varname (dict):    Dictionary variable-names on each mode diagnostic netcdf.
        dic_ncnames (dict):    Dictionary filenames of netcdf files per mode.
        modes (list):          Mode names (used as keys for dictionaries above).
        new_bins (list):       list or array with the bin limits.
        dic_files (dict):      Dictionary with info for output netcdf and base netcdf.
        newvarname (str):      New varname with diagnostic per bins.
        ftest (file, optional): Save test on opened file ftest. Defaults to None.
        save (bool, optional): Save netcdf (if not useful for test).
                               Defaults to True.
        test_info (dict)     : kind of test, area is a datarray of area_grid if needed by test.
                               Defaults to None.
                               Example: {kind='tendeny', area=datarray, years=['2000','2001']}

    Returns:
        None. It saves a netcdf.


    .. testcode::

        # EXAMPLE ================================================================
        # Case with: modes named by A and B  -> TWO MODES
        #            bins = [0.1, 5.0, 10.0] -> TWO BINS (three limits)
        #
        modes = ['A', 'B']
        bins  = [0.1, 5.0, 10.0]
        dic_fractions = {'A': [0.2, 0.5, 0.3], 'B': [0.1, 0.8, 0.1]}
        dic_varname = {'A': 'emidust', 'B': 'emidust' } 
        dic_ncnames = {'A': 'emidust_modeA_mymodel.nc', 'B': 'emidust_modeB_mymodel.nc' }
        dic_files   = {'basedir': '/home/myhome/mydirwithnc/',
                       'base_nc': 'my_canonical_netcdf.nc', 
                       # base_nc -> It could be any of the dic_ncnames.
                       'base_var': 'emidust',           
                       # base_bar -> It will be remove on save of binned file
                       'newf_nc':  'emidust_BINNED.nc'}
        newvarname  = 'emidust_bins'

        create_netcdf_emi_bins(dic_fractions, dic_varname, dic_ncnames, modes, new_bins,
                               dic_files,  newvarname, ftest=None, save=True)
  
        # (1) Read 'basedir + base_nc' -> copy from here the coords and dimensions of netcdf
        # (2) Add the bins dimensions to the given structure
        # (3) Open the files of dic_ncnames, for each file open the variable of dic_varname.
        # (4) With these files and the dic_fractions estimate the binned 2D file.
        # (5) As save=True it saved the netcdf with name newf_nc and variable name newvarname.
        #
        # For emissions or depositions it is a handy possibility to test consistency. 


        f_area = '/home/myhome/mydirwithnc/area_grid.nc'
        vararea = xr.open_dataset(f_area)['area']
        mytest = open('testing_results.txt', 'w')
        test_emi = {kind='tendency', area=vararea}

        create_netcdf_emi_bins(dic_fractions, dic_varname, dic_ncnames, modes, new_bins,
                               dic_files,  newvarname, ftest=mytest, save=True,
                               test_info=test_emi)

        mytest.close()


    .. testcode:

        # FULL EXAMPLE: Binning of emission fields of 4 dust modes of IPSL climate model.

        a_sigma = np.array([1.8, 2.0, 1.9, 2.00])
        a_mmd   = np.array([1.0, 2.5, 7.0, 22.0])
        a_mode  = ['D1','CI','SI','D4']
        new_bins = [0.0001, 0.2, 2.0, 3.6, 6.0, 12.0, 20.0, 30.0, 40.0, 100.0, 200.0, 300.0]

        dic_fractions = {}
        for sigma, mmd, mode in zip(a_sigma, a_mmd, a_mode): 
            fractions = ln.bin_fractions_lognormal(mmd, sigma, bins=new_bins)
            dic_fractions[mode] = fractions

        dic_varname = {'D1':'emidustD1','CI':'emidustCI', 'SI':'emidustSI','D4':'emidustD4'}

        dic_ncnames = {
                    'D1':'emidustD1_JKmon_IPSL-LMDZORINCAv6_4modes-DRE-global-ade_v1-r1i1p1f1_gr_20090101-20141231-accum.nc',
                    'CI':'emidustCI_JKmon_IPSL-LMDZORINCAv6_4modes-DRE-global-ade_v1-r1i1p1f1_gr_20090101-20141231-accum.nc',
                    'SI':'emidustSI_JKmon_IPSL-LMDZORINCAv6_4modes-DRE-global-ade_v1-r1i1p1f1_gr_20090101-20141231-accum.nc',
                    'D4':'emidustD4_JKmon_IPSL-LMDZORINCAv6_4modes-DRE-global-ade_v1-r1i1p1f1_gr_20090101-20141231-accum.nc'
                    }

        f_area = 'tests/area_grid.nc'
        vararea = xr.open_dataset(f_area)['area']

        test_emi = {'kind':'tendency', 'area':vararea, 'years':['2010','2011']}

        ftestemi = open('check_emission_binning.txt', 'w')

        dic_files = {'basedir': './tests',
                    'base_nc': dic_ncnames['D1'], 
                    'base_var': 'emidustD1',           
                    'newf_nc':  'emidust_BINNED.nc'}

        create_netcdf_2D_bins(dic_fractions, dic_varname, dic_ncnames, a_mode, new_bins,
                            dic_files,  'emidust_bin', ftest=ftestemi, save=True,
                            test_info=test_emi)

    
        # OUTPUT of check_emission_binning.txt ===================================


        ============== mode contributions ============================ 
        ---- 2010 ----------------------- 
            Contribution mode D1 :    97.58   |> accum =    97.58 
            Contribution mode CI :   715.22   |> accum =   812.80 
            Contribution mode SI :  5274.96   |> accum =  6087.76 
            Contribution mode D4 : 10673.08   |> accum = 16760.84 
            *Contribution ALL bins:                       16759.97 
        ---- 2011 ----------------------- 
            Contribution mode D1 :    88.92   |> accum =    88.92 
            Contribution mode CI :   651.81   |> accum =   740.73 
            Contribution mode SI :  4807.24   |> accum =  5547.97 
            Contribution mode D4 :  9726.73   |> accum = 15274.70 
            *Contribution ALL bins:                       15273.90 
        
        ==============  bins contributions ============================ 
        ---- 2010 ----------------------- 
            Contribution bin [   0.0,    0.2] :     0.40   |> accum =     0.40 
            Contribution bin [   0.2,    2.0] :   490.17   |> accum =   490.57 
            Contribution bin [   2.0,    3.6] :   946.51   |> accum =  1437.08 
            Contribution bin [   3.6,    6.0] :  1763.44   |> accum =  3200.52 
            Contribution bin [   6.0,   12.0] :  3858.86   |> accum =  7059.37 
            Contribution bin [  12.0,   20.0] :  3511.50   |> accum = 10570.87 
            Contribution bin [  20.0,   30.0] :  2635.21   |> accum = 13206.08 
            Contribution bin [  30.0,   40.0] :  1464.50   |> accum = 14670.57 
            Contribution bin [  40.0,  100.0] :  1935.78   |> accum = 16606.36 
            Contribution bin [ 100.0,  200.0] :   146.74   |> accum = 16753.10 
            Contribution bin [ 200.0,  300.0] :     6.87   |> accum = 16759.97 
            *Contribution ALL bins:                                   16759.97 
        ---- 2011 ----------------------- 
            Contribution bin [   0.0,    0.2] :     0.36   |> accum =     0.36 
            Contribution bin [   0.2,    2.0] :   446.71   |> accum =   447.07 
            Contribution bin [   2.0,    3.6] :   862.58   |> accum =  1309.66 
            Contribution bin [   3.6,    6.0] :  1607.08   |> accum =  2916.73 
            Contribution bin [   6.0,   12.0] :  3516.70   |> accum =  6433.44 
            Contribution bin [  12.0,   20.0] :  3200.15   |> accum =  9633.58 
            Contribution bin [  20.0,   30.0] :  2401.55   |> accum = 12035.13 
            Contribution bin [  30.0,   40.0] :  1334.64   |> accum = 13369.77 
            Contribution bin [  40.0,  100.0] :  1764.14   |> accum = 15133.91 
            Contribution bin [ 100.0,  200.0] :   133.73   |> accum = 15267.64 
            Contribution bin [ 200.0,  300.0] :     6.26   |> accum = 15273.90 
            *Contribution ALL bins:                                   15273.90 

    """


    base_dir = dic_files['basedir']   # where nc with modes are
    base_nc  = dic_files['base_nc']   # names of nc file as base file
    base_var = dic_files['base_var']  # varname in base_nc to replicate dims.
    newf_nc  = dic_files['newf_nc']   # new filename for netcdf binned
    
    bins_nc = xr.open_dataset(os.path.join(base_dir,base_nc)).copy(deep=True)
    bins_nc = bins_nc.assign_coords({'bins_up':   new_bins[1::]})
    bins_nc = bins_nc.assign_coords({'bins_down': new_bins[0:-1]})

    dic_nc_modes = {}
    dic_fullncnames = {mode:os.path.join(base_dir,dic_ncnames[mode]) for mode in modes}
    for mode in modes:
        dic_nc_modes[mode]=xr.open_dataset(dic_fullncnames[mode])[dic_varname[mode]]

    mode_shape = list(dic_nc_modes[modes[0]].shape)
    bins_shape = mode_shape+[int(len(new_bins)-1)]
    var_bin = np.zeros((*bins_shape,))  # addpt this with the info of mode_shape

    for kbin in range(len(new_bins)-1):
        for mode in modes:
            var_bin[:,:,:,kbin]=var_bin[:,:,:,kbin]+dic_nc_modes[mode]*dic_fractions[mode][kbin]

    newarray = xr.DataArray(var_bin, coords=[bins_nc['time'],bins_nc['lat'],
                                             bins_nc['lon'],bins_nc['bins_up']],
                            dims=['time','lat','lon','bins_up'])

    bins_nc = bins_nc.assign({newvarname:newarray})
    bins_nc = bins_nc.drop(base_var)

    if save==True:
        bins_nc.to_netcdf(os.path.join(base_dir,newf_nc)) 
    
    if test_info==None:
        return None
    elif test_info['kind']=='tendency':
        test_tendency(test_info, newarray, new_bins, modes,
                      base_dir, dic_fullncnames, dic_varname, ftest=ftest)
    elif test_info['kind']=='average':
        test_average()
    else:
        return
    
    return


def test_tendency(test_info, tendency_bins, new_bins, modes,
                   basedir, dic_ncnames, dic_varname, ftest=None):
    """
    This function compares a binned diagnostic :math:`X(t, \phi, \zeta; b)` with the set of files per mode :math:`\{X_{m}(t, \phi, \zeta)\}`, where b is the set of bins.

    It is assumed that test_info is of the form:
    test_info = { kind='tendency', area=datarray, years=list of years as str}

    For the calculation assumed a leap calendar and years with 365 days. Given that the test is 
    for consistency between modal and binned approaches last hypothesis is not important.

    Args:
        test_info (dict):         Dictionary for fractions per mode for new_bins.
        tendency_bins (datarray): Dictionary variable-names on each mode diagnostic netcdf.
        new_bins (list):          List or array with the bin limits.
        modes (list):             Mode names (used as keys for dictionaries above).
        dic_varname (dict):       Dictionary variable-names on each mode diagnostic netcdf.
        dic_ncnames (dict):       Dictionary filenames of netcdf files per mode.
        ftest (file, optional):    Save test on opened file ftest. Defaults to None.

    Returns:
        None. Print on screen or save to file.
    """

    vararea = test_info['area']
    lyears  = test_info['years']
 
    dic_checkbins = {}
    if ftest==None:
        print(' No save the double checking ')
        
    else:

        ftest.write(' \n ============== mode contributions ============================ \n')
        for year in lyears:
            tendency_test = tendency_bins.sum(['bins_up']).sel(time=year)
            check_bins = total_tendency(tendency_test, vararea)

            ftest.write('   ---- '+ year + ' ----------------------- \n')
            accum = 0
            for mode in modes:
                bins_nc_check = xr.open_dataset(dic_ncnames[mode])[dic_varname[mode]].sel(time=year)
                check_mode = total_tendency(bins_nc_check, vararea)
                accum = accum + check_mode
                ftest.write('     Contribution mode '+mode+' : %8.2f   |> accum = %8.2f \n' %(check_mode, accum))

            ftest.write(    '    *Contribution ALL bins:                       %8.2f \n' % (check_bins))

        ftest.write(' \n ==============  bins contributions ============================ \n')
        for year in lyears:
            tendency_test = tendency_bins.sum(['bins_up']).sel(time=year)
            check_bins = total_tendency(tendency_test, vararea)
            accum2= 0.0
            ftest.write('   ---- '+ year + ' ----------------------- \n')
            lcheck = []
            for bindown,binup in zip(new_bins[0:-1], new_bins[1::]):
                tendency_testbin = tendency_bins.sel(bins_up=str(binup)).sel(time=year)
                check_binup = total_tendency(tendency_testbin, vararea)
                lcheck.append(check_binup)
                accum2 = accum2 + check_binup
                ftest.write('     Contribution bin [%6.1f, %6.1f] : %8.2f   |> accum = %8.2f \n' %(bindown, binup, check_binup, accum2))
            dic_checkbins[year]=[a/accum2 for a in lcheck]
            ftest.write(    '    *Contribution ALL bins:                                    %8.2f \n' % (check_bins))
            
    return 

def _example_binning_emissions():
    """
    Testing example described in the documentation.

    It used the emissions files of 4 dust modes, and estimate
    the emissions per bin for a set of bins from 0.0001 micrometers to 300 micrometers. 
    """

    a_sigma = np.array([1.8, 2.0, 1.9, 2.00])
    a_mmd   = np.array([1.0, 2.5, 7.0, 22.0])
    a_mode  = ['D1','CI','SI','D4']
    new_bins = [0.0001, 0.2, 2.0, 3.6, 6.0, 12.0, 20.0, 30.0, 40.0, 100.0, 200.0, 300.0]

    dic_fractions = {}
    for sigma, mmd, mode in zip(a_sigma, a_mmd, a_mode): 
        fractions = ln.bin_fractions_lognormal(mmd, sigma, bins=new_bins)
        dic_fractions[mode] = fractions

    dic_varname = {'D1':'emidustD1','CI':'emidustCI', 'SI':'emidustSI','D4':'emidustD4'}

    dic_ncnames = {
                'D1':'emidustD1_JKmon_IPSL-LMDZORINCAv6_4modes-DRE-global-ade_v1-r1i1p1f1_gr_20090101-20141231-accum.nc',
                'CI':'emidustCI_JKmon_IPSL-LMDZORINCAv6_4modes-DRE-global-ade_v1-r1i1p1f1_gr_20090101-20141231-accum.nc',
                'SI':'emidustSI_JKmon_IPSL-LMDZORINCAv6_4modes-DRE-global-ade_v1-r1i1p1f1_gr_20090101-20141231-accum.nc',
                'D4':'emidustD4_JKmon_IPSL-LMDZORINCAv6_4modes-DRE-global-ade_v1-r1i1p1f1_gr_20090101-20141231-accum.nc'
                }

    f_area = 'tests/area_grid.nc'
    vararea = xr.open_dataset(f_area)['area']

    test_emi = {'kind':'tendency', 'area':vararea, 'years':['2010','2011']}

    ftestemi = open('check_emission_binning.txt', 'w')

    dic_files = {'basedir': './tests',
                'base_nc': dic_ncnames['D1'], 
                'base_var': 'emidustD1',           
                'newf_nc':  'emidust_BINNED.nc'}

    create_netcdf_2D_bins(dic_fractions, dic_varname, dic_ncnames, a_mode, new_bins,
                        dic_files,  'emidust_bin', ftest=ftestemi, save=True,
                        test_info=test_emi)

    return


if __name__ == "__main__":
    _example_binning_emissions()

