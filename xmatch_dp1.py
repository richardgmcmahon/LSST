"""xmatch with LSST DP1

https://sdm-schemas.lsst.io/dp1.html

Reference band - parameters measured on coadds of this band were used for
multi-band forced photometry


JWST image viewer:

https://s3.amazonaws.com/grizli-v2/ClusterTiles/Map/gds/jwst.html?coord=53.04105,-27.85449&zoom=5

http://grizli-cutout.herokuapp.com/

"""

# standard library functions
import os
import sys
import time

import glob

import configparser
import getpass
import inspect
from inspect import currentframe, getframeinfo
import logging
import traceback

import lsst_util as lu
# help(lu)

# 3rd party functions
import numpy as np
from scipy import stats

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import astropy
from astropy.table import Table, hstack, vstack
from astropy.coordinates import (search_around_sky, SkyCoord,
                                 match_coordinates_sky)
import astropy.units as u
from astropy.stats import mad_std
from astropy.io import ascii

from rgm_util import explore_flux_fluxerr
from rgm_util import plot_colour_magnitude
from rgm_util import xmatch_ls
from rgm_util import write_radec_csvfile
from rgm_util.mk_urls import *
from rgm_util import get_githash
from rgm_util import plot_radec
#from plot_radec import *

import astrolinks

# help(plot_radec)

import xmatch as xm
# help(xm)



def xmatch_tables(infile1=None,
                  infile2=None,
                  table1=None,
                  table2=None,
                  colnames_radec_table1=['RAJ2000', 'DECJ2000'],
                  colnames_radec_table2=['RAJ2000', 'DECJ2000'],
                  checkplots=True,
                  showplots=True,
                  show_table_in_browser=True,
                  verbose = False,
                  infostats=True,
                  selfmatch=False,
                  multimatch = False,
                  xmatch_seplimit_initial=2.0,
                  xmatch_seplimit_final=0.5,
                  xmatch_background_radius_limits=(4.0, 8.0)):
    """

    xmatch

    """

    from astropy.table import Table
    from astropy.table import hstack, vstack

    if infile1 is not None and table1 is None:
        table1 = Table.read(infile1)
        table1.meta['Filename'] = infile1
        logging.info(f'Table read in: {infile1}')

    infile1 = table1.meta['Filename']
    filename1 = os.path.basename(infile1)
    logging.info(f"Table1: {table1.meta['Filename']}")
    logging.info(f'Number of rows: {len(table1)}')
    logging.info(f'Number of columns: ' +
                 f'{len(table1.colnames)} {len(table1.columns)}')
    logging.info(f'Elapsed time(secs):  {time.time() - t0}\n')

    if infostats:
        table1.info(['attributes', 'stats'])
        logging.info(f'Elapsed time(secs): {time.time() - t0}\n')


    if infile2 is not None and table2 is None:
        table2 = Table.read(infile2)
        table2.meta['Filename'] = infile2
        logging.info(f'Table read in: {infile2}')

    infile2 = table2.meta['Filename']
    filename2 = os.path.basename(infile2)
    logging.info(f"Table2: {table2.meta['Filename']}")
    logging.info(f'Number of rows: {len(table2)}')
    logging.info(f'Number of columns: ' +
                 f'{len(table2.colnames)} {len(table2.columns)}')
    logging.info(f'Elapsed time(secs):  {time.time() - t0}\n')

    if infostats:
        table2.info(['attributes', 'stats'])
        logging.info(f'Elapsed time(secs): {time.time() - t0}\n')


    # window the xmatch secondary table/file2 based on
    # primary table1/infile1 RA, Dec limits to spped up the xmatch
    # This ignores 24hr wrap
    ra_limits1 = (np.min(table1[colnames_radec_table1[0]]),
                 np.max(table1[colnames_radec_table1[0]]))
    ra_range1 = ra_limits1[1] - ra_limits1[0]

    dec_limits1 = (np.min(table1[colnames_radec_table1[1]]),
                 np.max(table1[colnames_radec_table1[1]]))
    dec_range1 = dec_limits1[1] - dec_limits1[0]

    logging.info(f'RA range: {ra_limits1} {ra_range1}')
    logging.info(f'Dec range: {dec_limits1} {dec_range1}')

    logging.info(f'Elapsed time(secs): {time.time() - t0}\n')

    suptitle_prefix = filename1
    plotfile_prefix = filename1 + '_'
    logger.info('showplots: {showplots}')
    plot_radec(table=table1,
               colnames_radec=colnames_radec_table1,
               plotfile_prefix=plotfile_prefix,
               suptitle_prefix=suptitle_prefix,
               showplot=showplots)

    suptitle_prefix = filename2
    plotfile_prefix = filename2 + '_'
    plot_radec(table=table2,
               colnames_radec=colnames_radec_table2,
               plotfile_prefix=plotfile_prefix,
               suptitle_prefix=suptitle_prefix,
               showplot=showplots)


    tableWindow = True
    if tableWindow is True:
        table2_windowed = radec_window(
            table=table2,
            colnames_radec=colnames_radec_table2,
            ra_limits=ra_limits1,
            dec_limits=dec_limits1)

        nrows2 = len(table2_windowed)
        logging.info(
            f'Number of rows in windowed table2: {nrows2}')

        table2 = table2_windowed

        suptitle_prefix = filename2
        plotfile_prefix = filename2 + '_windowed_'
        plot_radec(table=table2,
               colnames_radec=colnames_radec_table2,
               plotfile_prefix=plotfile_prefix,
                   suptitle_prefix=suptitle_prefix,
                   showplot=showplots)


    nrows2 = len(table2)
    if nrows2 == 0:
        logging.warning(f'Exiting function since no rows in windowed table')
        results = None
        return results


    logging.info(f'xMatch primary reference file: {table1.meta["Filename"]}')
    logging.info(
        f'Number of rows per table1: {len(table1)} {len(table1)}')

    logging.info(f'xMatch secondary file: {table2.meta["Filename"]}')
    logging.info(
        f'Number of rows per table2: {len(table2)} {len(table2)}')
    logging.info(f'Initial Matching radius(arcsec): {xmatch_seplimit_initial}')
    logging.info(f'Elapsed time(sec): {time.time() - t0}')


    idx2, dr, dra, ddec = xm.xmatch_cat(
        table1=table1,
        table2=table2,
        colnames_radec1=colnames_radec_table1,
        colnames_radec2=colnames_radec_table2,
        seplimit=xmatch_seplimit_initial,
        multimatch=multimatch,
        xmatch_background_radius_limits=xmatch_background_radius_limits,
        verbose=verbose)

    logging.info(f'Crossmatch completed')
    logging.info(f'Elapsed time(secs): {time.time() - t0}\n')
    logging.info(f'{len(idx2)}, {len(dr)}')
    if not multimatch:
        logging.info(f'Number of matched sources: {len(table1)}, ' +
                     f'{len(table2)} {len(idx2)}')
    if multimatch:
        logging.info(f'Number of matched sources: {len(table1)}, ' +
                     f'{len(table2)}, {len(idx2[0])}, {len(idx2[1])}')

    logging.info(f'Maximum separation: {np.max(dr)}, {len(dr)}')
    logging.info(f'Delta RA range: {np.min(dra)}, {np.max(dra)}')
    logging.info(f'Delta Dec range: {np.min(ddec)} {np.max(ddec)}')

    seplimit = xmatch_seplimit_initial
    seplimit_test_scalefactors = [0.25, 0.50, 1.0, 2.0, 4.0]
    for seplimit_test_scalefactor in seplimit_test_scalefactors:
        seplimit_test = seplimit * seplimit_test_scalefactor
        itest = (dr < seplimit_test)
        print(seplimit_test, len(itest), len(dr[itest]), itest[0], itest[-1])
        if not multimatch:
            logging.info(f'Separation limit: {seplimit_test}, ' +
                     f'{len(table1[itest])}')
        if multimatch:
            logging.info(f'Separation limit: {seplimit_test}, ' +
                     f'{len(idx2[0][itest])}, {len(idx2[1][itest])}')

    itest = (dr <= seplimit)

    if checkplots:
        if not multimatch:
            ra1 = table1[itest][colnames_radec_table1[0]]
            dec1 = table1[itest][colnames_radec_table1[1]]
            ra2 = table2[idx2[itest]][colnames_radec_table2[0]]
            dec2 = table2[idx2[itest]][colnames_radec_table2[1]]
        if multimatch:
            ra1 =  table1[idx2[0][itest]][colnames_radec_table1[0]]
            dec1 = table1[idx2[0][itest]][colnames_radec_table1[1]]
            ra2 =  table2[idx2[1][itest]][colnames_radec_table2[0]]
            dec2 = table2[idx2[1][itest]][colnames_radec_table2[1]]

        plotfile_label = filename1 + '_x_' + filename2
        plotfile_prefix = filename1 + '_x_' + filename2 + '_'
        plotfile_suffix = None
        showplot = showplots
        xm.xmatch_checkplots(ra1=ra1, dec1=dec1,
                             ra2=ra2, dec2=dec2,
                             rmax=seplimit,
                             suptitle=plotfile_label,
                             title=plotfile_label,
                             plotfile_label=plotfile_label,
                             plotfile_prefix=plotfile_prefix,
                             plotfile_suffix=plotfile_suffix,
                             githash=False,
                             verbose=verbose,
                             showplot=showplot)

    seplimit_final = xmatch_seplimit_final
    itest = (dr < seplimit_final)

    if checkplots:
        logging.info(f'len(itest): {len(itest)}')
        logging.info(f'len(table1): {len(table1)}')
        logging.info(f'len(table2): {len(table2)}')
        logging.info(f'len(dr): {len(dr)}')
        logging.info(f'len(dra): {len(dra)}')
        logging.info(f'len(ddec): {len(ddec)}')
        if multimatch:
            logging.info(f'len(idx2[0]): {len(idx2[0])}')
            logging.info(f'len(idx2[1]): {len(idx2[1])}')
        if not multimatch:
            logging.info(f'len(idx2): {len(idx2)}')


        if not multimatch:
            ra1 = table1[itest][colnames_radec_table1[0]]
            dec1 = table1[itest][colnames_radec_table1[1]]
            ra2 = table2[idx2[itest]][colnames_radec_table2[0]]
            dec2 = table2[idx2[itest]][colnames_radec_table2[1]]
        if multimatch:
            ra1 =  table1[idx2[0][itest]][colnames_radec_table1[0]]
            dec1 = table1[idx2[0][itest]][colnames_radec_table1[1]]
            ra2 =  table2[idx2[1][itest]][colnames_radec_table2[0]]
            dec2 = table2[idx2[1][itest]][colnames_radec_table2[1]]



        plotfile_label = filename1
        plotfile_prefix = filename1

        xm.xmatch_checkplots(ra1=ra1, dec1=dec1,
                             ra2=ra2, dec2=dec2,
                             rmax=seplimit_final,
                             suptitle=None,
                             title=None,
                             plotfile_label=plotfile_label,
                             plotfile_prefix=plotfile_prefix,
                             plotfile_suffix=plotfile_prefix,
                             githash=False,
                             verbose=verbose,
                             showplot=showplot)

    # hstack the columns
    infostats = True
    saveresults = True
    if multimatch:
        result = hstack([table1[idx2[0]][itest],
                         table2[idx2[1]][itest]])

    if not multimatch:
        result = hstack([table1[itest], table2[itest]])


    result.add_columns([dr[itest],
                            dra[itest],
                            ddec[itest]],
                           names=['dr', 'dra', 'ddec'])


    if show_table_in_browser:
        result.show_in_browser(jsviewer=True)
        input('Enter any key to continue... ')

    result.write('tmp.html', format='ascii.html', overwrite=True)
    result.write('tmp2.html', format='jsviewer', overwrite=True)

    if infostats:
        result.info(['attributes', 'stats'])

    if saveresults:
        filename1 = os.path.basename(table1.meta['Filename'])
        filename2 = os.path.basename(table2.meta['Filename'])
        resultfile_prefix = filename1 + '_x_' + filename2
        outfile = resultfile_prefix + '.fits'
        result.write(outfile, overwrite=True)
        logging.info(f'Saving xmatch: {outfile}')


    # outer join
    print(len(dr), len(idx2))
    itest = (dr > seplimit)
    print(f'itest = (dr > seplimit) range: {np.min(itest)} {np.max(itest)}')
    dr[itest] = -999.9
    dra[itest] = -999.9
    ddec[itest] = -999.9

    nrows = len(table_lsst)
    redshift = np.empty(nrows, dtype=np.float32)
    #redshift[itest] = -999.9

    #result_outer = hstack([table1[itest],
    #                 table2[idx2[itest]]])




    return result

    # reverse cross-match
    table_desi = table1
    colnames_radec_table_desi = colnames_radec_table1

    table1 = table2
    table2 = table_desi

    colnames_radec1 = colnames_radec_table2
    colnames_radec2 = colnames_radec_table_desi
    logging.info(f'xMatch primary reference file: {table1.meta["Filename"]}')
    logging.info(f'xMatch secondary file: {table2.meta["Filename"]}')

    idx2, dr, dra, ddec = xm.xmatch_cat(
        table1=table1,
        table2=table2,
        colnames_radec1=colnames_radec1,
        colnames_radec2=colnames_radec2,
        seplimit=seplimit, verbose=verbose)

    logging.info(f'Elapsed time(secs): {time.time() - t0}')
    logging.info(f'Number of matched sources: ' +
                 f'{len(table1)} {len(table2)} {len(idx2)}')

    for seplimit_test_scalefactor in seplimit_test_scalefactors:
        seplimit_test = seplimit * seplimit_test_scalefactor
        itest = (dr < seplimit_test)
        logging.info(f'Separation limit: {seplimit_test} ' +
                     f'{len(table1[itest])}')

    itest = (dr < seplimit)
    if checkplots:

        ra1 = table1[itest][colnames_radec1[0]]
        dec1 = table1[itest][colnames_radec1[1]]
        ra2 = table2[idx2[itest]][colnames_radec2[0]]
        dec2 = table2[idx2[itest]][colnames_radec2[1]]

        plotfile_label = filename1 + '_x_' + filename2
        plotfile_prefix = filename1 + '_x_' + filename2
        xm.xmatch_checkplots(ra1=ra1, dec1=dec1,
                             ra2=ra2, dec2=dec2,
                             rmax=seplimit,
                             plotfile_label=plotfile_label,
                             plotfile_prefix=plotfile_prefix,
                             githash=False,
                             verbose=verbose)

        logging.info(f'Elapsed time(secs): time.time() - t0)')



    # hstack the columns
    result2 = hstack([table2[idx2[itest]],
                     table1[itest]])

    result2.add_columns([dr[itest],
                        dra[itest],
                        ddec[itest]],
                       names=['dr', 'dra', 'ddec'])

    if infostats:
        result2.info(['attributes', 'stats'])

    nrows1 = len(result1)
    ncolumns1 = len(result1.columns)
    logging.info(
        f'Number of matched sources and columns in result1: ' +
        f'{nrows1} {ncolumns1}')


    nrows2 = len(result2)
    ncolumns2 = len(result2.columns)
    logging.info(
        f'Number of matched sources and columns in result2: ' +
        f'{nrows2} {ncolumns2}')


    # now do a xmatch with all neighbours within a search radius
    table1 = table_desi
    table2 = table_secondary

    colnames_radec1 = colnames_radec_table_desi
    colnames_radec2 = colnames_radec_table_secondary

    idx2, dr, dra, ddec = xm.xmatch_cat(
        table1=table1,
        table2=table2,
        colnames_radec1=colnames_radec1,
        colnames_radec2=colnames_radec2,
        seplimit=seplimit,
        verbose=verbose)


    logging.info(f'Elapsed time(secs): {time.time() - t0}')
    logging.info(f'Number of matched sources: ' +
                 f'{len(table1)} {len(table2)} {len(idx2)}')


    results = [result1, result2]

    logging.info(
        f'Exiting {inspect.currentframe().f_code.co_name}')


    return results


def get_radec_limits(table=None,
                     colnames_radec=None):

    import numpy as np

    ra_limits = (np.min(table[colnames_radec[0]]),
                 np.max(table[colnames_radec[0]]))

    dec_limits = (np.min(table[colnames_radec[1]]),
                 np.max(table[colnames_radec[1]]))


    return ra_limits, dec_limits


def radec_window(table=None,
                 colnames_radec=['RAJ2000', 'DEJ2000'],
                 ra_limits = None, dec_limits = None,
                 ra_centre = None, dec_centre=None,
                 window_overlap= 0.01,
                 radius = None):
    """ All units are degrees """
    # add window edge around the RA, Dec limits with Cosine Dec correction
    # applied to RA limits.
    # np.cos(np.deg2rad(-60.00)) = 0.50000

    maxabsdec = np.max(np.abs(dec_limits))
    logging.info(f'Maximum abs(dec): {maxabsdec}')
    cosdec_correction = 1.0/(np.cos(np.deg2rad(maxabsdec)))
    logging.info(f'1/Cosine abs(dec): {cosdec_correction}')

    window_ra_limits = [ra_limits[0] - (window_overlap * cosdec_correction),
                        ra_limits[1] + (window_overlap * cosdec_correction)]
    window_dec_limits = [dec_limits[0] - window_overlap,
                         dec_limits[1] + window_overlap]

    print(f'window_ra_limits: {window_ra_limits}')
    print(f'window_dec_limits: {window_dec_limits}')

    ra_test = table[colnames_radec[0]]
    dec_test = table[colnames_radec[1]]
    itest = ((ra_test >= window_ra_limits[0]) &
             (ra_test <= window_ra_limits[1]) &
             (dec_test >= window_dec_limits[0]) &
             (dec_test <= window_dec_limits[1]))

    table_windowed = table[itest]

    return table_windowed




def explore_all(table=None, plot_title=None,
                plotfile_prefix=None,
                markersize=4.0):

    xrange = (-0.5, 1.0)
    yrange = (23.0, 14.0)

    lu.plot_cmodel_psf(table=table,
                       xrange=xrange,
                       yrange=yrange,
                       plot_title=plot_title,
                       plotfile_prefix=plotfile_prefix)


    lu.plot_cmodel_psf(table=table,
                       refExtendedness=True,
                       xrange=xrange,
                       yrange=yrange,
                       plot_title=plot_title,
                       plotfile_prefix=plotfile_prefix)


    return


def getargs():
    """

    """

    argfs = ""

    return args


def mk_mylogger(logfile_prefix=''):

    import sys
    import getpass
    import logging

    username = getpass.getuser()
    print('__name__:', __name__)

    logger = logging.getLogger(__name__)

    logfile_prefix = os.path.splitext(os.path.basename(__file__))[0]
    logfile = logfile_prefix + '_' + filename_timestamp + '.log'
    # set up file and screen logging
    # datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(name)-12s %(module)s %(funcName)s %(lineno)d %(levelname)-8s %(message)s',
        datefmt='%y-%m-%dT%H:%M:%S',
        filename=logfile,
        filemode='w')

    # define a Handler which writes INFO messages or higher to sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter(
        '%(asctime)s.%(msecs)03d %(name)-12s: %(module)s %(funcName)s %(lineno)s %(levelname)-8s %(message)s',
        datefmt='%y-%m-%dT%H:%M:%S')
    # tell the handler to use this format
    console.setFormatter(formatter)

    # add the handler to the root logger
    logging.getLogger().addHandler(console)

    logger.info('Username: ' + username)
    logger.info(__name__)
    logger.info(__file__)


    return logger




# do your work here
if __name__ == "__main__":

    #plt.figure(figsize=(10,5))

    # import TempleModels as tm
    # colors_dict = tm.get_colors_nested_dict()

    timestamp = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
    filename_timestamp = time.strftime('%Y%m%dT%H%M', time.gmtime())

    logger = mk_mylogger()

    t0 = time.time()

    debug = False
    verbose = False

    show_table_in_browser = False
    show_table_in_browser_jsviewer = False

    xmatch_selfmatch = False
    xmatch_multimatch = True
    xmatch_seplimit_initial = 5.0
    xmatch_seplimit_final = 1.0

    xmatch_checkplots = False
    xmatch_showplots = False

    infostats = False

    run_LRD = False
    run_JWST = False
    run_Milliquas = False
    run_DESI_AGNQSO_VAC = True

    # move these to explore_dp1
    run_cmodel_psf = False
    run_color_color = False
    run_color_mag = False

    zrange = (0.0, 6.0)
    zmin_list = 4.0

    config = configparser.ConfigParser()
    if run_Milliquas:
        configfile = 'Milliquas.cfg'
        logger.info('Read configfile: ' + configfile)

        config.read(configfile)

        sectionName = 'DEFAULT'
        inpath_milliquas  = config.get(sectionName, 'path')
        filename_milliquas  = config.get(sectionName, 'filename')
        infile_milliquas = inpath_milliquas + filename_milliquas


        infile_milliquas = 'DP1_Fields_Milliquas.fits'
        logging.info(f'Reading: {infile_milliquas}')
        milliquas = Table.read(infile_milliquas)
        milliquas.meta['Filename'] = infile_milliquas

        nrows = str(len(milliquas))
        logging.info('Number of rows:' + nrows)
        logging.info('Number of columns: ' +
                     str(len(milliquas.colnames)) + ' ' +
                     str(len(milliquas.columns)))
        if infostats:
            milliquas.info(['attributes', 'stats'])

        colnames_radec_milliquas = ['RAdeg', 'DEdeg']
        colname_redshift = 'z'

    if run_DESI_AGNQSO_VAC:
        infile_xmatch = 'DP1_Fields_DESI_DR1_AGNQSO_VAC.fits'
        logging.info(f'Reading: {infile_xmatch}')
        table_xmatch = Table.read(infile_xmatch)
        table_xmatch.meta['Filename'] = infile_xmatch

        # if infostats:
        table_xmatch.info()

        colnames_radec_xmatch = ['TARGET_RA', 'TARGET_DEC']
        colname_redshift = 'Z'

        table1 = table_xmatch
        colnames_radec_table1 = colnames_radec_xmatch


    if run_LRD:
        filename = 'LRD_Setton+2024.fits'
        infile_LRD = filename
        table = Table.read(infile_LRD)
        table.meta['Filename'] = filename
        table.info(['attributes', 'stats'])
        colnames_radec_lrd = ['RA', 'Dec']
        table_lrd = table

    if run_JWST:
        filename = 'jwst-sources.fits'
        infile = filename
        table = Table.read(infile)
        table.meta['Filename'] = filename
        table.info(['attributes', 'stats'])

        colnames_radec_xmatch = ['ra', 'dec']
        colname_redshift = 'zspec'

        zmin_list = 5.5
        zrange = (0.0, 7.0)

        table_xmatch = table


    inpath_lsst = './'
    # filename_lsst = 'dp1_Object_ecdfs.fits'
    # filename_lsst = 'DP1_ECDFS_Object.fits'
    # filename_lsst = 'objtab_dp1_EDFS.fits'
    # filename_lsst = 'objtab_dp1_LELF.fits'
    filename_lsst = 'DP1_ObjectThin.fits'



    infile_lsst = inpath_lsst + filename_lsst
    colnames_radec_lsst = ['coord_ra', 'coord_dec']

    logging.info(f'Reading infile1: {infile_lsst}')
    table_lsst = Table.read(infile_lsst)
    table_lsst.meta['Filename'] = infile_lsst
    nrows = str(len(table_lsst))
    logging.info('Number of rows: ' + nrows)
    logging.info('Number of columns: ' +
             str(len(table_lsst.colnames)) + ' ' +
             str(len(table_lsst.columns)))
    logging.info(f'Elapsed time(secs): {time.time() - t0}')

    if infostats:
        table_lsst.info(['attributes', 'stats'])
    logging.info(f'Elapsed time(secs): {time.time() - t0}\n')

    lu.count_refband(table=table_lsst)

    lu.explore_refExtendedness(table=table_lsst)

    logger.info(f"{table_lsst.meta['Filename']} {len(table_lsst)}")

    # explore_all(table=table_lsst)
    plot_title = table_lsst.meta['Filename']
    plotfile_prefix = table_lsst.meta['Filename'] + '_'

    wavebands = ['u', 'g', 'r', 'i', 'z', 'y']

    xrange = (-0.5, 1.0)
    yrange = (15.0, 25.0)


    """
    https://sdm-schemas.lsst.io/dp1.html#Object
     example Object flags

    {band}_blendedness_flag: https://sdm-schemas.lsst.io/dp1.html#Object.g_blendedness_flag
    {band_extendedness_flag: https://sdm-schemas.lsst.io/dp1.html#Object.g_extendedness_flag
    {band}_i_flag: https://sdm-schemas.lsst.io/dp1.html#Object.g_blendedness_flag


    """


    ra_limits_lsst, dec_limits_lsst = get_radec_limits(
        table=table_lsst, colnames_radec=colnames_radec_lsst)

    logging.info(f'ra_limits_lsst: {ra_limits_lsst}')
    logging.info(f'dec_limits_lsst: {dec_limits_lsst}')
    input('Enter any key to continue... ')

    if run_Milliquas:
        n_milliquas_all = len(milliquas)
        print(f'ra_limits_lsst: {ra_limits_lsst}')
        print(f'dec_limits_lsst: {dec_limits_lsst}')
        milliquas = radec_window(
            table=milliquas,
            colnames_radec=colnames_radec_milliquas,
            ra_limits=ra_limits_lsst, dec_limits=dec_limits_lsst)

        n_milliquas_windowed = len(milliquas)

        print(f'n_milliquas_all: {n_milliquas_all}')
        print(f'n_milliquas_windowed: {n_milliquas_windowed}')

        table1 = milliquas
        colnames_radec_table1 = colnames_radec_milliquas


    if run_LRD:
        table1 = table_lrd
        colnames_radec_table1 = colnames_radec_lrd

    if run_JWST:
        table1 = table_xmatch
        colnames_radec_table1 = colnames_radec_xmatch
        xmatch_AN = True
        xmatch_multimatch = True

    print('Now do the xmatch')
    input('Enter any key to continue... ')

    table2 = table_lsst
    colnames_radec_table2 = colnames_radec_lsst

    ra_limits1, dec_limits1= get_radec_limits(
        table=table_lsst,
        colnames_radec=colnames_radec_lsst)

    logging.info(f'RA limits: {ra_limits1[0]}, {ra_limits1[1]}')
    ra_range1 = ra_limits1[1] - ra_limits1[0]
    logging.info(f'RA range: {ra_range1}\n')

    logging.info(f'Dec limits: {dec_limits1[0]}, {dec_limits1[1]}')
    dec_range1 = dec_limits1[1] - dec_limits1[0]
    logging.info(f'DEC range: {dec_range1}\n')


    logging.info(f"xmatch table1: {table1.meta['Filename']}; {len(table1)}")
    logging.info(f'{colnames_radec_table1}')
    logging.info(f"xmatch table2: {table2.meta['Filename']}; {len(table2)}")
    logging.info(f'{colnames_radec_table2}')
    logging.info(f'xmatch selfmatch: {xmatch_selfmatch}')
    logging.info(f'xmatch multimatch: {xmatch_multimatch}')
    logging.info(f'xmatch initial radial seplimit (") {xmatch_seplimit_initial}')
    logging.info(f'xmatch final radial seplimit (") {xmatch_seplimit_final}')
    input('Enter any key to continue... ')

    table_xmatch = xmatch_tables(
        table1=table1,
        table2=table2,
        colnames_radec_table1=colnames_radec_table1,
        colnames_radec_table2=colnames_radec_table2,
        xmatch_seplimit_initial=xmatch_seplimit_initial,
        xmatch_seplimit_final=xmatch_seplimit_final,
        selfmatch=xmatch_selfmatch,
        multimatch=xmatch_multimatch,
        xmatch_background_radius_limits=(5.0, 10.0),
        checkplots=xmatch_checkplots,
        showplots=xmatch_showplots,
        show_table_in_browser=show_table_in_browser,
        verbose=verbose,
        infostats=False)

    table_xmatch.info(['attributes', 'stats'])

    filename1 = os.path.basename(table1.meta['Filename'])
    filename2 = os.path.basename(table2.meta['Filename'])
    plotfile_prefix = filename1 + '_x_' + filename2 + '_'
    plot_title = filename1 + '_x_' + filename2 + '\n' \
        + f'r_sep <= {xmatch_seplimit_final:.2f}"'

    table = table_xmatch

    logging.info(f'{type(table)}')
    logging.info(f'{len(table)}')


    # plot_title = 'dp1_Object_ecdfs.fits x Milliquas'
    # plot redshift histogram

    xcolname = colname_redshift
    xdata= table[xcolname]

    itest_zgte = (xdata >= zmin_list)
    print(f'Number with z>= {zmin_list}; {len(table[itest_zgte])}')
    if show_table_in_browser_jsviewer:
        table[itest_zgte].show_in_browser(jsviewer=True)
    input('Enter any key to continue... ')

    for itest, test_z in enumerate(itest_zgte):
        if test_z:
            print()
            print(test_z, xdata[itest])
            print(test_z, table[itest])
            print(test_z, table['coord_ra'][itest])
            print(test_z, table['coord_dec'][itest])
            url_ls = astrolinks.mk_link_LS(
                RA=table['coord_ra'][itest],
                Dec=table['coord_dec'][itest], DR='dr10')
            print(f'Legacy Survey viewer link: {url_ls}')

    print(f'Number with z>= {zmin_list}: {len(table[itest_zgte])}')
    logger.info('')
    input('Enter any key to continue... ')

    # sys.exit()

    plotfile_prefix = os.path.basename(table1.meta['Filename']) + \
        '_xm_' + \
        os.path.basename(table2.meta['Filename'])

    xcolname = colname_redshift
    lu.plot_hist_redshift(table=table,
                          colname_redshift=colname_redshift,
                          bins=30, zrange=(0.0, 6.0),
                          plot_title=plot_title,
                          plotfile_prefix=plotfile_prefix)

    if run_DESI_AGNQSO_VAC:
        logging.info('\n')
        lu.desi_plot_hist_redshift(table=table,
                                   colname_redshift=colname_redshift,
                                   plot_title=plot_title,
                                   plotfile_prefix=plotfile_prefix)

    markersize=6.0
    plotfile_prefix = os.path.basename(table1.meta['Filename']) + \
        '_xm_' + \
        os.path.basename(table2.meta['Filename'])

    xrange = (-0.5, 1.0)
    yrange = (15.0, 25.0)

    logging.info('\n')
    lu.plot_cmodel_psf(table=table,
                       xrange=xrange,
                       yrange=yrange,
                       plot_title=plot_title,
                       markersize=markersize,
                       plotfile_prefix=plotfile_prefix)

    logging.info('\n')
    lu.plot_cmodel_psf(table=table,
                       refExtendedness=True,
                       xrange=xrange,
                       yrange=yrange,
                       plot_title=plot_title,
                       markersize=markersize,
                       plotfile_prefix=plotfile_prefix)






    lu.plot_redshift_mag(table=table,
                         magtype='psfMag',
                         markersize=markersize,
                         colname_redshift=colname_redshift,
                         zrange=zrange,
                         plot_title=plot_title)

    plt.figure(figsize=(10,6))
    lu.plot_redshift_color(table=table,
                           magtype='psfMag',
                           markersize=markersize,
                           colname_redshift=colname_redshift,
                           zrange=zrange,
                           plot_title=plot_title)
