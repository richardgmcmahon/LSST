"""
https://sdm-schemas.lsst.io/dp1.html


"""

import os
import sys
import time
import logging
import configparser

import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, join_skycoord, join

global DEBUG

import TempleModels as tm
#from TempleModels import * as tm
# help(TempleModels)

# help(TempleModels.plot_ugr)

from lsst.utils.plotting import (get_multiband_plot_colors,
                                 get_multiband_plot_symbols,
                                 get_multiband_plot_linestyles)

filter_colors = get_multiband_plot_colors()
filter_symbols = get_multiband_plot_symbols()
filter_linestyles = get_multiband_plot_linestyles()

dp1_FieldRadius = 1.0 # degree

field_filters = ['u', 'g', 'r', 'i', 'z', 'y']
bands = field_filters

magrange = (15.0, 26.0)

logger = logging.getLogger(__name__)


import argparse

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='lsst utilities')
    parser.add_argument('-d' , '--debug',
                        action='store_true',
                        help='Enable debug mode')
    args = parser.parse_args()
    return args


def get_dp1_fieldnames():
    """Retrieve LSST DP1 field names from configuration file.
    
    Reads field names for LSST Data Preview 1 (DP1) from the
    lsst_dp1.cfg configuration file.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    list of str
        List of DP1 field names (e.g., ['47Tuc', 'LELF', 'FRSG', 'ECDFS', 'EDFS']).
    
    Notes
    -----
    The field names correspond to specific sky regions observed
    in the LSST DP1 dataset. The configuration file must exist
    in the current working directory.
    """

    logger.info('\n')

    config = configparser.ConfigParser()

    configfile = 'lsst_dp1.cfg'
    sectionName = 'DP1'

    print(f'Read configfile: {configfile} {sectionName}')
    logger.info(f'Read configfile: {configfile} {sectionName}')

    config.read(configfile)
    fieldnames = config.get(sectionName, 'fieldnames')
    fieldnames = [item.strip() for item in fieldnames.split(',')]
    print(f'fieldnames: {fieldnames}')

    """
    # Convert the string to a Python list
    my_list = [item.strip() for item in my_list_str.split(',')]

    print(my_list)
    # Output: ['apple', 'banana', 'cherry', 'date']
    """

    # fieldnames = ['47Tuc', 'LELF', 'FRSG', 'ECDFS', 'EDFS', 'LGLF', 'SeaNeb']


    return fieldnames


def read_Milliquas_DP1():
    """Read Milliquas quasar catalog matched to LSST DP1 fields.
    
    Loads the Milliquas quasar catalog that has been pre-matched to
    LSST Data Preview 1 field regions.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    ra_MQ : astropy.table.Column
        Right ascension in degrees (J2000).
    dec_MQ : astropy.table.Column
        Declination in degrees (J2000).
    redshift_MQ : astropy.table.Column
        Redshift values.
    
    Notes
    -----
    Reads from file 'DP1_Fields_Milliquas.fits' in the current directory.
    The Milliquas catalog is a comprehensive compilation of QSOs and AGN.
    """

    logger.info('\n')

    infile = 'DP1_Fields_Milliquas.fits'
    milliquas = Table.read(infile)
    milliquas.info(['attributes', 'stats'])
    ra_MQ = milliquas['RADeg']
    dec_MQ = milliquas['DEDeg']
    redshift_MQ = milliquas['z']

    return ra_MQ, dec_MQ, redshift_MQ


def xmatch_Milliquas(table_lsst=None,
                     table_MQ=None,
                     searchRadius=0.5):
    """Cross-match LSST catalog with Milliquas quasar catalog.
    
    Performs spatial cross-matching between an LSST photometric catalog
    and the Milliquas quasar compilation.
    
    Parameters
    ----------
    table_lsst : astropy.table.Table, optional
        LSST source catalog with RA/Dec columns.
    table_MQ : astropy.table.Table, optional
        Milliquas catalog table with RA/Dec and redshift columns.
    searchRadius : float, optional
        Matching radius in arcseconds. Default is 0.5.
    
    Returns
    -------
    index_lsst : int
        Index placeholder for LSST matches.
    index_MQ : int
        Index placeholder for Milliquas matches.
    
    Notes
    -----
    This is a placeholder function that needs full implementation.
    """
    logger.info('\n')

    index_lsst = 1
    index_MQ = 1

    return index_lsst, index_MQ



def count_refband(table=None):
    """Count objects by reference band in LSST catalog.
    
    Analyzes the distribution of objects across different reference
    photometric bands in an LSST Object catalog.
    
    Parameters
    ----------
    table : astropy.table.Table
        LSST Object catalog with 'refBand' column.
    
    Returns
    -------
    None
        Prints band counts to stdout.
    
    Notes
    -----
    The reference band is the photometric band used for
    multi-band forced photometry in the LSST pipeline.
    Typically one of: u, g, r, i, z, y.
    """
    logger.info('')

    nrows = len(table)
    unique, unique_counts = \
        np.unique(table['refBand'], return_counts=True)
    total_counts = np.sum(unique_counts)

    print(unique)
    print(unique_counts)
    print(f'refBand counts: {total_counts} {nrows}')

    for index, row in enumerate(unique):
        print(index, row, unique_counts[index])
    print(f'Total number of source: {len(table)}')
    print()

    return


def plot_redshift_mag(table=None,
                      magtype='psfMag',
                      zrange=(0.0, 6.0),
                      magrange=(15.0, 25.0),
                      markercolor=None,
                      markersize=4.0,
                      colname_redshift='z',
                      figsize=(8.0, 6.0),
                      plot_title=None,
                      plotfile_prefix='',
                      plotpath=None,
                      showplot=True,
                      savefig=True):
    """

    """
    # plot redshift v mag hubble diagram

    logger.info('\n')

    bands = ['u', 'g','r', 'i', 'z', 'y']
    magtype = 'psfMag'

    if plot_title is None:
        plot_title = table.meta['Filename']

    xcolname = colname_redshift
    for band in bands[0:6]:

        plt.figure(figsize=figsize)

        xdata_save = table[xcolname]

        ycolname = band + '_' + magtype
        ydata_save = table[ycolname]

        ndata_finite = np.count_nonzero(
            ~np.isnan(table[band + '_extendedness']))
        table[band + '_extendedness'].info(['attributes', 'stats'])
        print(f'Number of finite values: {ndata_finite}')

        label = str(ndata_finite)
        plt.plot(np.nan, np.nan, color='None', label=label)

        itest = (table[band + '_extendedness'] == 1)
        xdata = xdata_save[itest]
        ydata = ydata_save[itest]
        label = str(len(xdata)) + ' Extended'
        plt.plot(xdata, ydata, '.r',
                 markersize=markersize,
                 label=label)

        xdata= table[xcolname]
        ydata = table[ycolname]

        itest = (table[band + '_extendedness'] == 0)
        xdata = xdata_save[itest]
        ydata = ydata_save[itest]
        label = str(len(xdata)) + ' Pointlike'
        plt.plot(xdata, ydata, '.b', label=label)

        plt.xlim(zrange)
        plt.ylim(magrange)

        plt.legend()
        plt.xlabel(xcolname)
        plt.ylabel(ycolname)
        plt.title(plot_title)
        plt.gca().invert_yaxis()

        plt.show()

    return


def table_flag_clean(table=None,
                     bandlist=None,
                     band_flaglist=None,
                     flaglist=None):
    """Filter LSST catalog by removing flagged sources.
    
    Removes sources with quality flags set in specified bands,
    returning a cleaned catalog.
    
    Parameters
    ----------
    table : astropy.table.Table
        LSST source catalog with flag columns.
    bandlist : list of str, optional
        List of photometric bands to check (e.g., ['u', 'g', 'r']).
        Default is ['u', 'g', 'r', 'i', 'z', 'y'].
    band_flaglist : list of str, optional
        List of band-specific flag column suffixes to check
        (e.g., ['extendedness_flag', 'blendedness_flag']).
    flaglist : list of str, optional
        List of non-band-specific flag columns to check.
    
    Returns
    -------
    astropy.table.Table
        Filtered table with flagged sources removed (flag == 0 retained).
    
    Notes
    -----
    Flags are combined across all specified bands and flag types.
    Only sources with flag value of 0 (good) are retained.
    """
    logger.info('\n')

    if bandlist is None:
        bandlist = ['u', 'g','r', 'i', 'z', 'y']

    for band in bandlist[0:6]:

        if band_flaglist is not None:
            for band_flag in band_flaglist:
                itest = (table[band + '_' + band_flag] == 0)
                table = table[itest]

    return table


def table_flag_info(table=None,
                    bandlist=None,
                    band_flaglist=None,
                    flaglist=None):
    """Report statistics on LSST catalog quality flags.
    
    Prints information about the number of sources passing
    quality flag criteria in each photometric band.
    
    Parameters
    ----------
    table : astropy.table.Table
        LSST source catalog with flag columns.
    bandlist : list of str, optional
        List of photometric bands to check (e.g., ['u', 'g', 'r']).
        Default is ['u', 'g', 'r', 'i', 'z', 'y'].
    band_flaglist : list of str, optional
        List of band-specific flag column suffixes to report
        (e.g., ['extendedness_flag', 'blendedness_flag']).
    flaglist : list of str, optional
        List of non-band-specific flag columns to report.
    
    Returns
    -------
    None
        Prints flag statistics to stdout and log.
    
    Notes
    -----
    Reports counts of sources with flag == 0 (good quality)
    for each combination of band and flag type.
    """

    logger.info('\n')

    if bandlist is None:
        bandlist = ['u', 'g','r', 'i', 'z', 'y']

    table_save = table


    for band in bandlist[0:6]:

        table = table_save
        print()
        logging.info(f'Total number of sources {len(table)}')
        if band_flaglist is not None:
            for band_flag in band_flaglist:
                itest = (table[band + '_' + band_flag] == 0)
                ntrue = itest.sum()
                logging.info(
                    f'Sources with {band}_{band_flag} == 0: ' +
                    f'{ntrue}')
                table = table[itest]

    return


def explore_extendedness(table=None,
                         band='ref',
                         figsize=(12, 6),
                         plot_suptitle=None,
                         plotfile_prefix=None,
                         showplot=True,
                         plotpath=None):
    """

    Explore refExtendedness and refSizeExtendedness

    https://sdm-schemas.lsst.io/dp1.html

    refExtendedness
    https://sdm-schemas.lsst.io/dp1.html#Object.refExtendedness

    refSizeExtendedness
    https://sdm-schemas.lsst.io/dp1.html#Object.refSizeExtendedness
    Moments-based measure of whether an object is point-like (0) or extended (1). Reference band.

    """
    logger.info('\n')

    plt.close(plt.gcf())

    debug = True

    filename = os.path.basename(table.meta['Filename'])

    colname_list = ['refBand', 'refExtendedness', 'refSizeExtendedness']
    table[colname_list].info(['attributes', 'stats'])
    logger.info('\n')

    table['refBand'].info(['attributes', 'stats'])
    data = table['refBand']
    has_mask = hasattr(data, 'mask')
    n_masked = 0
    if has_mask:
        imasked = np.ma.is_masked(data)
        n_masked = len(data[imasked])
        print(f'n_masked {n_masked}')
        n_masked = np.sum(data.mask)
    print(f'n_masked {n_masked}')
    logger.info('\n')

    table['refExtendedness'].info(['attributes', 'stats'])
    logger.info('\n')

    table['refSizeExtendedness'].info(['attributes', 'stats'])
    logger.info('\n')

    fig, axes = plt.subplots(1, 2, figsize=(12,6))

    colnames_subplots = ['refExtendedness', 'refSizeExtendedness']
    for iplot, colname in enumerate(colnames_subplots):

        table[colname].info(['attributes', 'stats'])

        xdata = table[colname]
        if hasattr(xdata, "mask"):
            n_bad = np.count_nonzero(xdata.mask)
        else:
            try:
                n_bad = np.count_nonzero(np.isinf(xdata) | np.isnan(xdata))
            except Exception:
                n_bad = 0
        ndata = len(xdata)
        ndata_bad = n_bad
        ndata_good = ndata - ndata_bad

        legend_title = 'All/ Good/ Masked'
        print(f'ndata all, good, bad  = {ndata} {ndata_good} {ndata_bad}')
        label = f'{ndata}/ {ndata_good}/ {ndata_bad}'
        plot_title = filename

        axes[iplot].hist(xdata, bins=100, label=label)
        axes[iplot].set_yscale('log')

        axes[iplot].legend(title=legend_title, loc='upper center')
        axes[iplot].set_xlabel(colname)
        axes[iplot].set_ylabel('Number per bin')
        axes[iplot].set_title(plot_title)

    if plot_suptitle is not None:
        fig.suptitle(plot_suptitle)

    timestamp = time.strftime('%Y-%m-%dT%H:%M', time.gmtime())

    footnote1 = timestamp + ': ' + \
            os.path.basename(__file__)

    plt.figtext(0.01, 0.01, footnote1,
                    ha="left", fontsize=8,
                    bbox={"facecolor":'none', 'edgecolor':'none',
                          "alpha":0.5, "pad":2})

    footnote2 = os.path.basename(__file__)
    plt.figtext(0.99, 0.50, footnote2,
                    rotation='vertical',
                    ha="right", fontsize=8,
                    bbox={"facecolor":'none', 'edgecolor':'none',
                          "alpha":0.5, "pad":2})

    logger.info(f'plotfile_prefix: {plotfile_prefix}')
    if plotfile_prefix is None:
        plotfile_prefix = filename

    plotfile = (plotfile_prefix +
                '_hist_' + colnames_subplots[0] +
                '_' + colnames_subplots[1] + '.png')

    print(f'Saving plotfile: {plotfile}')
    plt.savefig(plotfile)
    logger.info('\n')
    if showplot:
        plt.show()

    print(f'bands {bands}')
    print('Loop through the bands')
    if debug:
        input('Enter any key to continue... ')

    for band in bands:
        fig, axes = plt.subplots(1, 2, figsize=(12,6))
        colnames_subplots = [band + '_extendedness',
                             band + '_sizeExtendedness']

        for iplot, colname in enumerate(colnames_subplots):

            table[colname].info(['attributes', 'stats'])

            xdata = table[colname]
            ndata = len(xdata)

            if hasattr(xdata, "mask"):
                n_bad = np.count_nonzero(xdata.mask)
            else:
                try:
                    n_bad = np.count_nonzero(np.isinf(xdata) | np.isnan(xdata))
                except Exception:
                    n_bad = 0

            ndata = len(xdata)
            ndata_bad = n_bad
            ndata_good = ndata - ndata_bad

            legend_title = 'All/ Good/ Masked'
            print(f'ndata all, good, bad  = {ndata} {ndata_good} {ndata_bad}')
            label = f'{ndata}/ {ndata_good}/ {ndata_bad}'
            plot_title = table.meta['Filename']

            if ndata_good > 0:
                axes[iplot].hist(xdata, bins=100, label=label)
                axes[iplot].set_yscale('log')

                axes[iplot].legend(title=legend_title, loc='upper center')
                axes[iplot].set_xlabel(colname)
                axes[iplot].set_ylabel('Number per bin')
                axes[iplot].set_title(plot_title)

        if plot_suptitle is not None:
            fig.suptitle(plot_suptitle)

        logger.info('\n')
        logger.info(f'plotfile_prefix: {plotfile_prefix}')

        if plotfile_prefix is None:
            plotfile_prefix = filename

        plotfile = (plotfile_prefix +
                '_hist_' + colnames_subplots[0] +
                '_' + colnames_subplots[1] + '.png')

        print(f'Saving plotfile: {plotfile}')
        plt.savefig(plotfile)
        logger.info('\n')

        if showplot:
            plt.show()

    logger.info('\n')

    colname_list = ['refBand', 'refExtendedness']
    for colname in colname_list:
        logging.info(f'{colname} {len(table)}')
        unique, unique_counts = \
            np.unique(table[colname], return_counts=True)

        print(unique)
        print(unique_counts)
        for index, row in enumerate(unique):
            print(index, row, unique_counts[index])
        logger.info('\n')

    logger.info('\n')



    # see https://docs.astropy.org/en/latest/_modules/astropy/utils/data_info.html#DataInfo
    data = table['refExtendedness']
    if hasattr(data, "mask"):
        n_bad = np.count_nonzero(data.mask)
    else:
        try:
            n_bad = np.count_nonzero(np.isinf(data) | np.isnan(data))
        except Exception:
            n_bad = 0
    ndata_bad = n_bad
    print(f'ndata_bad = {ndata_bad}')

    ndata_inf = np.count_nonzero(np.isinf(table['refExtendedness']))
    print(f'ndata_inf = {ndata_inf}')
    ndata_nan = np.count_nonzero(np.isnan(table['refExtendedness']))
    print(f'ndata_nan = {ndata_nan}')

    ndata_finite = np.count_nonzero(
            ~np.isnan(table['refExtendedness']))

    ndata = len(table)
    ndata_notfinite = ndata - ndata_finite
    print(f'Number of finite refExtendedness values : {ndata_finite}')
    print(f'Number of not finite refExtendedness values : {ndata_notfinite}')
    logger.info('\n')

    return


def plot_cmodel_psf(table=None,
                    band=None,
                    band_flaglist=None,
                    refExtendedness=False,
                    markersize=4.0,
                    xrange=None,
                    yrange=None,
                    multiplot=False,
                    figsize=(7.0, 7.0),
                    plot_title=None,
                    plotfile_prefix='',
                    savefig=True,
                    showplot=True,
                    debug=False):
    """

    plot cmodelMag - psfMag


    """
    # global logger
    logging.info('Starting ')

    debug = True


    plt.figure(figsize=figsize)

    if plot_title is None:
        plot_title = table.meta['Filename']

    if band is None:
        bands = ['u', 'g','r', 'i', 'z', 'y']
    if band is not None:
        bands = band

    magtype1 = 'psfMag'
    magtype2 = 'cModelMag'

    logger.info('\n')
    print('Next cycle through bands')
    if debug:
        input('Enter any key to continue... ')


    table_save = table
    for band in bands[0:6]:

        logger.info(f'band: {band}')

        plt.figure(figsize=figsize)

        table = table_save

        logger.info(f'band_flaglist: {band_flaglist}')
        if band_flaglist is not None:
            for band_flag in band_flaglist:
                itest = (table[band + '_' + band_flag] == 0)
                table = table[itest]

        print()
        logging.info('')
        legend_title = band + '_extendedness'

        ndata = len(table)
        # print(table[band + '_extendedness'])

        colname_extendedness = band + '_extendedness'
        table[colname_extendedness].info(['attributes', 'stats'])
        logger.info('\n')

        data = table[colname_extendedness]
        if hasattr(data, "mask"):
           n_bad = np.count_nonzero(data.mask)
        else:
            try:
               n_bad = np.count_nonzero(np.isinf(data) | np.isnan(data))
            except Exception:
               n_bad = 0
        ndata_bad = n_bad
        print(f'ndata_bad = {ndata_bad}')

        ndata_inf = np.count_nonzero(np.isinf(table[colname_extendedness]))
        print(f'ndata_inf = {ndata_inf}')
        ndata_nan = np.count_nonzero(np.isnan(table[colname_extendedness]))
        print(f'ndata_nan = {ndata_nan}')


        ndata_finite = np.count_nonzero(
            ~np.isnan(table[colname_extendedness]))

        print(f'Number of finite {band} ExtendedNess values : {ndata_finite}')
        print()

        ndata_finite = np.count_nonzero(
            ~np.isnan(table[band + '_psfMag']))
        table[band + '_psfMag'].info(['attributes', 'stats'])
        print(f'Number of finite {band}_psfMag values : {ndata_finite}')


        ndata_finite = np.count_nonzero(
            ~np.isnan(table[band + '_cModelMag']))
        table[band + '_cModelMag'].info(['attributes', 'stats'])
        print(f'Number of finite {band}_cModelMag values : {ndata_finite}')
        print()

        unique, unique_counts = \
            np.unique(table[band + '_extendedness'], return_counts=True)
        print(band + '_extendedness', len(table))
        print(unique)
        print(unique_counts)
        for index, row in enumerate(unique):
            print(index, row, unique_counts[index])
        print()

        data = table[band + '_extendedness']
        if hasattr(data, "mask"):
            n_bad = np.count_nonzero(data.mask)
        else:
            try:
                n_bad = np.count_nonzero(np.isinf(data) | np.isnan(data))
            except Exception:
                n_bad = 0
        ndata_bad = n_bad
        logger.info(f'ndata_bad = {ndata_bad}')

        label = (str(ndata) + '/ ' + str(ndata_finite) +
                 '/ ' + str(ndata_bad))
        plt.plot(np.nan, np.nan, color='None', label=label)

        xlabel = band + '_' + magtype1 + ' - ' + \
            band + '_' + magtype2
        xdata_save = table[band + '_' + magtype1] - \
            table[band + '_' + magtype2]

        ndata_finite = np.count_nonzero(~np.isnan(xdata_save))
        print(f'Number of finite {band}_{magtype1} - {band}_{magtype2} values : {ndata_finite}')

        logger.info('\n')

        ycolname = band + '_' + magtype1
        ydata_save = table[ycolname]

        itest = (table[colname_extendedness] != 1) & \
            (table[colname_extendedness] != 0)

        extendedness = table[colname_extendedness]
        print(f'extendness range: ' +
              f'{np.min(extendedness[itest])} ' +
              f'{np.max(extendedness[itest])} ' +
              f'{np.median(extendedness[itest])}')

        xdata = xdata_save[itest]
        ydata = ydata_save[itest]
        ndata = len(xdata)
        label = str(ndata)
        plt.plot(xdata, ydata, 'ko',
                 markersize=markersize/2.0,
                 markerfacecolor='none',
                 markeredgecolor='black',
                 label=label)


        # extended sources
        if not refExtendedness:
            itest = (table[band + '_extendedness'] == 1)
        if refExtendedness:
            itest = (table['refExtendedness'] == 1)

        xdata = xdata_save[itest]
        ydata = ydata_save[itest]
        label = str(len(xdata)) + ' Extended'
        plt.plot(xdata, ydata, '.r',
                 markersize=markersize,
                 markeredgecolor='none',
                 label=label)

        if not refExtendedness:
            itest = (table[band + '_extendedness'] == 0)
        if refExtendedness:
            itest = (table['refExtendedness'] == 0)
            legend_title = 'refExtendedness'

        xdata = xdata_save[itest]
        ydata = ydata_save[itest]
        label = str(len(xdata)) + ' Pointlike'
        plt.plot(xdata, ydata, '.b',
                 markersize=markersize,
                 label=label,
                 markeredgecolor='none')

        plt.xlim(-0.5, 2.5)
        if xrange is not None:
            plt.xlim(xrange)

        plt.ylim(16.0, 28.0)
        if yrange is not None:
            plt.ylim(yrange)

        plt.gca().invert_yaxis()

        plt.legend(title=legend_title, loc='upper right')
        plt.grid()

        plt.xlabel(xlabel)
        plt.ylabel(ycolname)

        plt.title(plot_title)
        # plt.suptitle(plot_suptitle)

        timestamp = time.strftime('%Y-%m-%dT%H:%M', time.gmtime())
        footnote2 = os.path.basename(__file__) + ': ' + timestamp
        plt.figtext(0.99, 0.50, footnote2,
                rotation='vertical',
                ha="right", fontsize=8,
                bbox={"facecolor":'none', 'edgecolor':'none',
                      "alpha":0.5, "pad":2})


        if savefig:
            plotfile_suffix = band + '_' + magtype1 + '-' + \
            band + '_' + magtype2
            if refExtendedness:
                plotfile_suffix = plotfile_suffix + '_refExtendedness'
            if not refExtendedness:
                plotfile_suffix = plotfile_suffix + '_' \
                    + band + '_extendedness'
            print(f'plotfile_prefix: {plotfile_prefix}')
            print(f'plotfile_suffix: {plotfile_suffix}')
            plotfile = plotfile_prefix + plotfile_suffix + '.png'
            # logger.info(f'Saving plotfile: {plotfile}')
            print(f'Saving plotfile: {plotfile}')
            plt.savefig(plotfile)


        logging.info('')
        if showplot:
            plt.show()


    # 2up plot
    #plt.figure(figsize=(12.0, 6.0))
    #plt.subplot(121)

    #xdata = table[band + '_' + magtype1] - \
    #    table[band + '_' + magtype2]



    return


def plot_redshift_color(table=None,
                        bands = ['u', 'g','r', 'i', 'z', 'y'],
                         magtype='psfMag',
                         zrange=(0.0, 6.0),
                         markersize=4.0,
                         colname_redshift='z',
                         figsize=(8.0, 6.0),
                         plot_title=None,
                         plotfile_prefix='',
                         savefig=True):
    # redshift versus color

    logger.info('\n')

    if plot_title is None:
        plot_title = table.meta['Filename']

    xcolname = colname_redshift

    for iband, band in enumerate(bands[0:5]):

        plt.figure(figsize=figsize)

        xdata_save = table[xcolname]

        band1 = bands[iband]
        mag1 = band1 + '_' + magtype

        band2 = bands[iband + 1]
        mag2 = band2 + '_' + magtype

        ycolname = mag1 + ' - ' + mag2
        ydata_save = table[mag1] - table[mag2]
        ydata = ydata_save

        ndata_finite = np.count_nonzero(~np.isnan(ydata))
        print(f'Number of finite y-axis values: {ndata_finite}')

        label = str(ndata_finite)
        plt.plot(np.nan, np.nan, color='None', label=label)


        itest = (table[band2 + '_extendedness'] == 1)
        xdata = xdata_save[itest]
        ydata = ydata_save[itest]
        label = str(len(xdata)) + ' Extended'
        plt.plot(xdata, ydata, '.r', label=label,
                 markersize=markersize)

        itest = (table[band2 + '_extendedness'] == 0)
        xdata = xdata_save[itest]
        ydata = ydata_save[itest]
        label = str(len(xdata)) + ' Pointlike'
        plt.plot(xdata, ydata, '.b', label=label,
                 markersize=markersize)

        plt.xlim(zrange)
        plt.ylim(-1.5, 5.0)

        plt.title(plot_title)
        plt.legend()
        plt.xlabel(xcolname)
        plt.ylabel('Reshift ' + ycolname)

        logging.info('')
        plt.show()

    return


def plot_color_mag(
        table=None,
        survey='LSST',
        band_ref=None,
        band1=None,
        band2=None,
        magtype=None,
        multiplot=True,
        colname_band_ref=None,
        colname_band1=None,
        colname_band2=None,
        label_band1='',
        label_band2='',
        label_band1_prefix='',
        label_band2_prefix='',
        label_band1_suffix='',
        label_band2_suffix='',
        xlabel=None,
        ylabel=None,
        xrange=(-2.5, 5.0),
        yrange=(25.0, 15.0),
        markercolor='black',
        title=None,
        suptitle=None,
        outpath='',
        plotpath='',
        filename_ls='',
        filename_prefix='',
        plotfile_prefix='',
        plotfile_suffix='',
        subplot2=True,
        showplots=True,
        overplot_TempleModel=False):
    """Plot color versus magnitude figure

    default is to column names for xlabel and ylabels


    """

    logger.info('\n')

    filename = os.path.basename(table1.meta['Filename'])
    print('plotfile_prefix:', plotfile_prefix)
    print('filename', filename)
    plot_title = filename

    print('Number of rows:', len(table))
    ndata_band1_finite = np.count_nonzero(np.isfinite(table[colname_band1]))
    ndata_band2_finite = np.count_nonzero(np.isfinite(table[colname_band2]))
    print('Finite band1 data:', ndata_band1_finite)
    print('Finite band2 data:', ndata_band2_finite)

    # 2up plot
    plt.figure(figsize=(12.0, 6.0))
    plt.subplot(121)

    xdata = table[colname_band1] - table[colname_band2]
    xlabel = colname_band1 + ' - ' + colname_band2

    ydata = table[colname_band1]
    ylabel = colname_band1


    if colname_band_ref is not None:
        ydata = table[colname_band_ref]
        ylabel = colname_band_ref

    if colname_band_ref is None:
        ydata = table[colname_band2]
        ylabel = colname_band2

    # define colname_ref for plot label usage; a bit messy
    if colname_band_ref is None:
        colname_band_ref = colname_band2

    ndata = len(xdata)
    nxdata_finite = np.count_nonzero(np.isfinite(xdata))
    nydata_finite = np.count_nonzero(np.isfinite(ydata))
    ndata_finite = np.count_nonzero(
            (np.isfinite(xdata) & np.isfinite(ydata)))

    label = (str(ndata) + '/' +
             str(nxdata_finite) + '/' +
             str(nydata_finite) + '/' +
             str(ndata_finite))
    print('All data:', ndata)
    print('Finite xdata:', nxdata_finite)
    print('Finite ydata:', nydata_finite)
    print('Finite data:', ndata_finite)
    print('label:', label)

    plt.plot(xdata, ydata, '.', ms=1.0,
             color=markercolor,
             label=label)
    plt.legend(loc='upper left')

    plt.gca().invert_yaxis()

    plt.suptitle(suptitle)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    plt.subplot(122)
    plt.plot(xdata, ydata, '.', ms=1.0,
             color=markercolor,
             label=label)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='upper left')

    plt.gca().invert_yaxis()
    plt.xlim((-2.5, 5.0))
    plt.ylim((25.0, 15.0))
    plt.title(title)


    # Legacy Survey
    if format == 'LS':
        plt.subplot(122)
        # itest = (table['EXT'] == 0) & (table['type'] == 'PSF')
        itest = (table['type'] == 'PSF')
        xdata = xdata[itest]
        ydata = ydata[itest]

        ndata = len(xdata)
        ndata_finite = np.count_nonzero(
            (np.isfinite(xdata) & np.isfinite(ydata)))
        label = str(ndata) + '/' + str(ndata_finite)
        plt.plot(xdata, ydata, '.r', ms=1.0,
             label=label)

        plt.legend(loc='upper left')

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.gca().invert_yaxis()
        plt.xlim(xrange)
        plt.ylim(yrange)

        plt.suptitle(plot_title)


    plotfile = outpath + (plotfile_prefix +
                         colname_band1 + '-' + colname_band2 +
                         '_v_' + colname_band1 + '.png')

    timestamp = time.strftime('%Y-%m-%dT%H:%M', time.gmtime())
    footnote1 = timestamp + ': ' + \
        os.path.basename(__file__) + ': ' + \
        plotfile
    plt.figtext(0.01, 0.01, footnote1,
                ha="left", fontsize=8,
                bbox={"facecolor":'none', 'edgecolor':'none',
                      "alpha":0.5, "pad":2})

    footnote2 = os.path.basename(__file__)
    plt.figtext(0.99, 0.50, footnote2,
                rotation='vertical',
                ha="right", fontsize=8,
                bbox={"facecolor":'none', 'edgecolor':'none',
                      "alpha":0.5, "pad":2})


    #    lineInfo()
    logger.info('\n')
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    if showplots:
        plt.show()

    return





def plot_color_magnitude2(table=None,
                           colname_redshift=None,
                           mag1='g',
                           mag2='r',
                           mag1_prefix='LS8_',
                           mag2_prefix='LS8_',
                           xlimits=(-1.0, 5.0),
                           ylimits=(25.0, 15.0),
                           title=None,
                           suptitle=None,
                           plotdir='',
                           plotfile_prefix='',
                           plotfile_suffix='',
                           showplots=True,
                           overplot=False):
    """
    2up plot m1 - m2 v m1 and m1 - m2 v m2

    """

    logger.info('\n')

    # color magnitude
    xdata = table[band1 + '_' + magtype] - table[band2 + '_' + magtype]
    ydata = data[band2 + '_' + magtype]

    xlabel = mag1_prefix + mag1 + ' - ' + \
        mag2_prefix  + mag2 + ' [AB] [DECam LS DR8]'
    ylabel = mag2_prefix + mag2 + ' [AB]'

    redshift = data[colname_redshift]
    print(len(data), len(redshift))
    itest = (redshift > 0.001)
    xdata = xdata[itest]
    ydata = ydata[itest]

    ndata = len(xdata)
    label = str(ndata)
    plt.plot(xdata, ydata,
             '.', ms=1.0,
             label=label)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if title is not None:
        plt.title(title, fontsize='small')

    if suptitle is not None:
        plt.suptitle(suptitle, fontsize='small')

    if overplot is False:
        plt.xlim(xlimits)
        plt.ylim(ylimits)

    ReadTempleModels.plot_color_magnitude(
        colors_dict=templemodels_dict,
        colormags=[mag1, mag2],
        mag=mag2,
        overplot=True)

    plt.legend()

    plotfile =  plotfile_prefix + \
        'ModelColors_Temple+2021_DECam_' + mag1 + '_' + mag2 + \
        '_' + mag2 + plotfile_suffix + '.png'
    if plotdir is not None:
        plotfile = plotdir + plotfile

    logger.info('\n')
    print('Saving plotfile:', plotfile)
    plt.savefig(plotfile)

    if showplots:
        plt.show()

    return


def plot_color_color(
        table=None,
        survey=None,
        colname_band_ref=None,
        colname_xband1=None,
        colname_xband2=None,
        colname_yband1=None,
        colname_yband2=None,
        colormagsx = None,
        colormagsy = None,
        xlabel='',
        ylabel='',
        xrange=(-2.5, 5.0),
        yrange=(-2.5, 5.0),
        markercolor='black',
        markersize=1.0,
        title=None,
        suptitle=None,
        outpath='',
        plotpath='',
        filename_ls='',
        filename_prefix='',
        plotfile_prefix='',
        plotfile_suffix='',
        subplot2=True,
        showplot=True,
        overplot_TempleModels=False,
        templemodels_dict=None):
    """Plot color versus magnitude figure

    default is to column names for xlabel and ylabels


    """
    import time
    from datetime import datetime

    logging.info('')

    if markercolor is None:
        markercolor = 'black'

    timestamp = time.strftime('%Y-%m-%dT%H:%M', time.gmtime())
    timestamp = datetime.now().strftime('%y%m%d_%H%M%S')


    filename = os.path.basename(table1.meta['Filename'])
    print('plotfile_prefix:', plotfile_prefix)
    print('filename', filename)
    plot_title = title
    if plot_title is None:
        plot_title = filename

    print('Number of rows:', len(table))
    ndata_xband1_finite = np.count_nonzero(np.isfinite(table[colname_xband1]))
    ndata_xband2_finite = np.count_nonzero(np.isfinite(table[colname_xband2]))
    print('Finite xband1 data:', ndata_xband1_finite)
    print('Finite xband2 data:', ndata_xband2_finite)

    # 2up plot
    plt.figure(figsize=(12.0, 6.0))
    plt.subplot(121)

    xdata = table[colname_xband1] - table[colname_xband2]
    xlabel = colname_xband1 + ' - ' + colname_xband2

    ydata = table[colname_yband1] - table[colname_yband2]
    ylabel = colname_yband1 + ' - ' + colname_yband2

    ndata = len(xdata)
    nxdata_finite = np.count_nonzero(np.isfinite(xdata))
    nydata_finite = np.count_nonzero(np.isfinite(ydata))
    ndata_finite = np.count_nonzero(
            (np.isfinite(xdata) & np.isfinite(ydata)))

    label = (str(ndata) + '/' +
             str(nxdata_finite) + '/' +
             str(nydata_finite) + '/' +
             str(ndata_finite))
    print('All data:', ndata)
    print('Finite xdata:', nxdata_finite)
    print('Finite ydata:', nydata_finite)
    print('Finite data:', ndata_finite)
    print('label:', label)

    plt.plot(xdata, ydata, '.',
             ms=markersize,
             color=markercolor,
             label=label)

    plt.legend(loc='upper left')

    plt.title(plot_title)
    plt.suptitle(suptitle)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


    plt.subplot(122)
    plt.plot(xdata, ydata, '.',
             ms=1.0,
             color=markercolor,
             label=label)

    plt.legend(loc='upper left')
    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.title(plot_title)
    plt.suptitle(suptitle)


    # Legacy Survey
    if format == 'LS':
        plt.subplot(122)
        # itest = (table['EXT'] == 0) & (table['type'] == 'PSF')
        itest = (table['type'] == 'PSF')
        xdata = xdata[itest]
        ydata = ydata[itest]

        ndata = len(xdata)
        ndata_finite = np.count_nonzero(
            (np.isfinite(xdata) & np.isfinite(ydata)))
        label = str(ndata) + '/' + str(ndata_finite)

        plt.plot(xdata, ydata, '.',
                 ms=markersize,
                 color=markercolor,
                 label=label)

        plt.legend(loc='upper left')

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.xlim(xrange)
        plt.ylim(yrange)

        plt.title(plot_title)
        plt.suptitle(plot_title)


    plotfile = outpath + (plotfile_prefix +
                         colname_xband1 + '-' + colname_xband2 +
                         '_v_' + colname_yband1 + '.png')

    timestamp = time.strftime('%Y-%m-%dT%H:%M', time.gmtime())
    footnote1 = timestamp + ': ' + \
        os.path.basename(__file__) + ': ' + \
        plotfile
    plt.figtext(0.01, 0.01, footnote1,
                ha="left", fontsize=8,
                bbox={"facecolor":'none', 'edgecolor':'none',
                      "alpha":0.5, "pad":2})

    footnote2 = os.path.basename(__file__) + ': ' + timestamp
    plt.figtext(0.99, 0.50, footnote2,
                rotation='vertical',
                ha="right", fontsize=8,
                bbox={"facecolor":'none', 'edgecolor':'none',
                      "alpha":0.5, "pad":2})


    #    lineInfo()
    logger.info('\n')
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    if showplot:
        plt.show()

    return


def plot_magnitude_distribution(table=None,
                                bands=None,
                                magbins=None,
                                suptitle=None,
                                grid=False,
                                extended_magtype=None):
    """
    From Tutorial

    301.4. Extended Chandra Deep Field South (ECDFS)

    https://dp1.lsst.io/tutorials/notebook/301/notebook-301-4.html

    Plot histograms of object magnitudes, separating the samples into stars
    and galaxies.

    Use the refExtendedness flag to distinguish them: treat objects with
    refExtendedness == 1 as likely galaxies (extended) and those with
    refExtendedness == 0 as likely stars (point sources).

    Use cModelMag mags for galaxies and psfMag for stars.

    """
    logger.info('\n')

    objtab = table

    field_filters = ['u', 'g', 'r', 'i', 'z', 'y']

    ptsource = (objtab['refExtendedness'] == 0)

    mag_bins = np.arange(16.0, 28.0, 0.2)

    nrows, ncols = 2, 3
    fig, ax = plt.subplots(
        nrows, ncols, figsize=(12, 8),
        sharey=True)
    ax = ax.flatten()
    plt.subplots_adjust(hspace=0, wspace=0)

    print()
    for iband, band in enumerate(field_filters):
        mag_psf = objtab[f"{band}_psfMag"]
        mag_cmodel = objtab[f"{band}_cModelMag"]

        valid_mags = ((15 < mag_psf) & (mag_psf < 30) &
                      (15 < mag_cmodel) & (mag_cmodel < 30))

        print(f"Number of objects in {band}-band: {np.sum(valid_mags)}")
        if np.sum(valid_mags) > 0:
            row = iband // ncols
            col = iband % ncols
            xdata = mag_psf[valid_mags & ptsource]
            label = band + ':' + str(len(xdata))
            print(f'iband: {iband}')
            n_psf, bins_psf, patches_psf = \
                ax[iband].hist(xdata,
                         bins=mag_bins,
                         histtype='step', linewidth=1,
                         color=filter_colors[band], label=label)
            for patch in patches_psf:
                patch.set_linestyle(filter_linestyles[band])

            ax[iband].set_xlim(mag_bins.min(), mag_bins.max())
            ax[iband].set_yscale('log')
            ax[iband].set_xlabel('Magnitude')
            if col == 0:
                ax[iband].set_ylabel('Number of objects')
            # ax[iband].set_title('Point sources (PSF mags)')

            ax[iband].minorticks_on()
            if grid:
                ax[iband].grid()

            xdata = mag_cmodel[valid_mags & ~ptsource]
            label = band + ':' + str(len(xdata))
            n_cmodel, bins_cmodel, patches_cmodel = \
                ax[iband].hist(xdata,
                         bins=mag_bins,
                         histtype='step',
                         linewidth=2,
                         color=filter_colors[band],
                         label=label,)
            for patch in patches_cmodel:
                patch.set_linestyle(filter_linestyles[band])

            xdata = mag_psf[valid_mags & ~ptsource]
            label = band + ':' + str(len(xdata))
            n_cmodel, bins_cmodel, patches_cmodel = \
                ax[iband].hist(xdata,
                               bins=mag_bins,
                               histtype='step',
                               linewidth=1,
                               color=filter_colors[band],
                               label=label,)
            for patch in patches_cmodel:
                    patch.set_linestyle(filter_linestyles[band])

            ax[iband].legend(loc='upper left', ncols=1)

    if suptitle is not None:
        plt.suptitle(suptitle)

    logger.info('\n')
    plt.show()

    return

    print()
    if extended_magtype != 'psf':
        fig, (ax1, ax2) = plt.subplots(
            1, 2,figsize=(10, 6),
            sharey=True)

    if extended_magtype == 'psf':
        fig, (ax1, ax2, ax3) = plt.subplots(
            1, 3, figsize=(14, 6),
            sharey=True)

    plt.subplots_adjust(hspace=0, wspace=0)

    for band in field_filters:
        mag_psf = objtab[f"{band}_psfMag"]
        mag_cmodel = objtab[f"{band}_cModelMag"]

        valid_mags = ((15 < mag_psf) & (mag_psf < 30) &
                      (15 < mag_cmodel) & (mag_cmodel < 30))

        print(f"Number of objects in {band}-band: {np.sum(valid_mags)}")
        if np.sum(valid_mags) > 0:
            xdata = mag_psf[valid_mags & ptsource]
            label = band + ':' + str(len(xdata))
            n_psf, bins_psf, patches_psf = \
                ax1.hist(xdata,
                         bins=mag_bins,
                         histtype='step', linewidth=2,
                         color=filter_colors[band], label=label)
            for patch in patches_psf:
                patch.set_linestyle(filter_linestyles[band])

            xdata = mag_cmodel[valid_mags & ~ptsource]
            label = band + ':' + str(len(xdata))
            n_cmodel, bins_cmodel, patches_cmodel = \
                ax2.hist(xdata,
                         bins=mag_bins,
                         histtype='step',
                         linewidth=2,
                         color=filter_colors[band],
                         label=label,)
            for patch in patches_cmodel:
                patch.set_linestyle(filter_linestyles[band])


            if extended_magtype == 'psf':
                xdata = mag_psf[valid_mags & ~ptsource]
                label = band + ':' + str(len(xdata))
                n_cmodel, bins_cmodel, patches_cmodel = \
                    ax3.hist(xdata,
                             bins=mag_bins,
                             histtype='step',
                             linewidth=2,
                             color=filter_colors[band],
                             label=label,)
                for patch in patches_cmodel:
                    patch.set_linestyle(filter_linestyles[band])

    ax1.set_xlim(mag_bins.min(), mag_bins.max())
    ax1.set_yscale('log')
    ax1.set_xlabel('Magnitude (psf)')
    ax1.set_ylabel('Number of objects')
    ax1.set_title('Point sources (PSF mags)')
    ax1.legend(loc='upper left', ncols=1)
    ax1.minorticks_on()
    if grid:
        ax1.grid()


    ax2.set_xlim(mag_bins.min(), mag_bins.max())
    ax2.set_xlabel('Magnitude (cModel)')
    ax2.set_title('Extended sources (cModel mags)')
    ax2.legend(loc='upper left', ncols=1)
    ax2.minorticks_on()
    if grid:
        ax2.grid()


    if extended_magtype == 'psf':
        ax3.set_xlim(mag_bins.min(), mag_bins.max())
        ax3.set_xlabel('Magnitude (psf)')
        ax3.set_title('Extended sources (psf mags)')
        ax3.legend(loc='upper left', ncols=1)
        ax3.minorticks_on()
        if grid:
            ax3.grid()

    if suptitle is not None:
        plt.suptitle(suptitle)

    logger.info('\n')
    plt.show()

    return



def plot_cmd_ccd(table=None,
                 bands=None,
                 suptitle=None,
                 cmap=None,
                 hexbins='log'):
    """
    Color magnitude and Color Color diagram

    Plot vs. CMDs and vs. CCDs. Create separate plots for stars and galaxies
    using PSF magnitudes for stars and cModel magnitudes for galaxies.

    bins = None produces linear scaling


    """
    import sys

    logger.info('\n')


    debug = locals().get('debug', False)
    if debug:
        try:
            help(plt.hexbin)
        except Exception as e:
            print(f"Warning: {e}")

    bins = hexbins
    if hexbins is not None:
        bins = hexbins

    print(f'hexbin bining: {bins}')
    objtab = table

    ptsource = (objtab['refExtendedness'] == 0)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        2, 2, figsize=(7, 6), height_ratios=[1.0, 1.5])
    plt.subplots_adjust(hspace=0, wspace=0)

    grmin, grmax = -1.0, 2.5
    rimin, rimax = -1.5, 3.0
    magmin, magmax = 16.0, 27.0

    bands = ['g', 'r', 'i']

    xdata = (objtab[ptsource][bands[0] + '_psfMag'] -
             objtab[ptsource][bands[1] + '_psfMag'])
    ydata = (objtab[ptsource][bands[1] + '_psfMag'] -
             objtab[ptsource][bands[2] + '_psfMag'])
    label = str(len(xdata))
    ax1.hexbin(xdata, ydata,
        gridsize=(150, 150),
        extent=(grmin, grmax, rimin, rimax),
        bins=bins,
        cmap='Grays', label=label)
    ax1.set_title('Point sources (PSF mags)')
    ax1.set_xlim(grmin, grmax)
    ax1.set_ylim(rimin, rimax)
    ax1.set_ylabel(r'$(r-i)$')
    ax1.set_xticklabels([])
    ax1.minorticks_on()
    ax1.annotate(label, xy=(0.1,0.9),
                 xycoords='axes fraction',
                 fontsize=12)

    cmap = 'Grays'
    cmap = plt.cm.Reds
    xdata = objtab[~ptsource]['g_cModelMag']-objtab[~ptsource]['r_cModelMag']
    ydata = objtab[~ptsource]['r_cModelMag']-objtab[~ptsource]['i_cModelMag']
    label = str(len(xdata))
    ax2.hexbin(xdata, ydata,
        gridsize=(150, 150),
        extent=(grmin, grmax, rimin, rimax), bins=bins,
        cmap=cmap, label=label)
    ax2.set_title('Extended sources (cModel mags)')
    ax2.set_xlim(grmin, grmax)
    ax2.set_ylim(rimin, rimax)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.minorticks_on()
    ax2.annotate(label, xy=(0.1,0.9),
                 xycoords='axes fraction',
                 fontsize=12)

    xdata = objtab[ptsource]['g_psfMag']-objtab[ptsource]['r_psfMag']
    ydata = objtab[ptsource]['r_psfMag']
    label = str(len(xdata))
    ax3.hexbin(xdata, ydata,
               gridsize=(150, 200),
        extent=(grmin, grmax, magmin, magmax), bins=bins,
        cmap='Grays')

    ax3.set_xlim(grmin, grmax)
    ax3.set_ylim(magmax, magmin)
    ax3.set_xlabel(r'$(g-r)$')
    ax3.set_ylabel(r'$r$ magnitude')
    ax3.minorticks_on()
    ax3.annotate(label, xy=(0.1,0.9),
                 xycoords='axes fraction',
                 fontsize=12)

    xdata = objtab[~ptsource]['g_cModelMag']-objtab[~ptsource]['r_cModelMag']
    ydata = objtab[~ptsource]['r_cModelMag']
    label = str(len(xdata))
    print(f'xdata range: {np.min(xdata)} {np.max(xdata)}')
    ax4.hexbin(xdata, ydata,
        gridsize=(150, 200),
        extent=(grmin, grmax, magmin, magmax), bins=bins,
        cmap='Reds')
    ax4.set_xlim(grmin, grmax)
    ax4.set_ylim(magmax, magmin)
    ax4.set_xlabel(r'$(g-r)$')
    ax4.set_yticklabels([])
    ax4.minorticks_on()
    ax4.annotate(label, xy=(0.1,0.9),
                 xycoords='axes fraction',
                 fontsize=12)


    if suptitle is not None:
        plt.suptitle(suptitle)

    logger.info('\n')
    plt.show()

    return


def plot_image_quality():
    """Plot image quality metrics for LSST observations.
    
    Placeholder function for generating plots of image quality
    metrics such as PSF FWHM, sky background, limiting magnitude,
    and airmass as a function of time or other parameters.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
    
    Notes
    -----
    This function needs to be implemented with actual plotting
    functionality for LSST image quality metrics.
    """

    return



def explore_column(data=None,
                   table=None,
                   colname=None):
    """ see also table_stats

    a column is specified as data or a table; could add pandas option

    """

    if data is not None:
        print(type(data))

    if table is not None:
        print(type(table))

    logger = logging
    logger.info('\n')

    if table is not None:
        table[colname].info(['attributes', 'stats'])
        data = table[colname]
        print(f'colname: {colname}')
        logger.info(f'colname: {colname}')

    print(type(data))
    is_masked = np.ma.is_masked(data)
    print(f'Data is_masked: {is_masked}')

    if is_masked:
        data_masked = data[data.mask]
        data_unmasked = data[~data.mask]
        print('Number of rows; all, unmasked, masked: ' +
          f'{len(data)} {len(data_unmasked)} {len(data_masked)}')

    print(f'data[0, -1]: {data[0]} {data[-1]}')

    ndata_finite = np.count_nonzero(np.isfinite(data))
    print(f'Number of finite values: {ndata_finite}')

    itest_isnan  = np.isnan(data)
    print('Number of NaN values:', len(data[itest_isnan]))
    itest_isinf  = np.isinf(data)
    print('Number of INF values:', len(data[itest_isinf]))
    itest_isfinite  = np.isfinite(data)
    print('Number of finite values:', len(data[itest_isfinite]))
    print()

    imasked = np.ma.is_masked(data)
    print('imasked:', imasked)
    if imasked:
        itest_is_mask = np.ma.is_mask(data.mask)
        itest_is_masked = data.mask

        print('Numpy Masked Array min, max, median:',
              np.ma.min(data), np.ma.max(data),
              np.ma.median(data))

        print('Masked values:', np.ma.count_masked(data),
              np.ma.count(data),
              len(data[itest_is_mask]),
              len(data[itest_is_masked]))

    ndata_finite = np.count_nonzero(np.isfinite(data))
    print(f'All: {len(data)}')
    print(f'Finite and not NaN: {ndata_finite}')

    print(f'min, max, median: {np.min(data)} ' +
          f'{np.max(data)} ' +
          f'{np.median(data)} ' +
          f'{len(data)}')
    print(f'nan min, max, median: {np.nanmin(data)} ' +
          f'{np.nanmax(data)} ' +
          f'{np.nanmedian(data)} ' +
          f'{len(data)}')
    print(f'masked min, max, median: {np.ma.min(data)} ' +
          f'{np.ma.max(data)} ' +
          f'{np.ma.median(data)} ' +
          f'{len(data)}')

    unique, unique_counts = \
        np.unique(data, return_counts=True)
    total_counts = np.sum(unique_counts)
    print(f'Maximum unique counts: {np.max(unique_counts)}')


    return


def  plot_hist_psfFlux_S_N(table=None,
                           fluxtype='psfFlux',
                           fieldname=None,
                           flags=None):
    """

    https://sdm-schemas.lsst.io/dp1.html#Object

    """

    logging.info('logging')
    logger.info('logger')


    fluxtype_options = ['psfFlux', 'free_psfFlux', 'cModelFlux',
                        'free_cModelFlux']

    plt.figure(figsize=(8, 6))
    linewidth = 2

    band = 'refBand'
    print(f'fluxtype: {fluxtype}')
    logger.info(f'fluxtype: {fluxtype}')
    colname = band + '_' + fluxtype + '_S_N'
    stats_column(table=table, colname=colname)

    xdata = table[colname]
    ndata_finite = np.count_nonzero(np.isfinite(xdata))

    ndata_S_N_gte_5 = np.count_nonzero(xdata >= 5.0)
    ndata_S_N_lt_0 = np.count_nonzero(xdata < 0.0)
    label = ('ref ' + str(len(xdata)) +
             ' ;' + str(ndata_finite) +
             ' ;' + str(ndata_S_N_gte_5) +
             ' ;' + str(ndata_S_N_lt_0))
    xrange = (-5.0, 25.0)
    n, bins, patches = plt.hist(xdata, bins=30,
                                histtype='step',
                                linewidth=linewidth,
                                range=xrange, label=label)


    for band in bands:

        colname = band + '_' + fluxtype + '_S_N'
        stats_column(table=table, colname=colname)
        xdata = table[colname]

        ndata_finite = np.count_nonzero(np.isfinite(xdata))
        ndata_S_N_gte_5 = np.count_nonzero(xdata >= 5.0)
        ndata_S_N_lt_0 = np.count_nonzero(xdata < 0.0)
        xrange = (-5.0, 25.0)

        label = (band + '   ' + str(len(xdata)) +
                 ' ;' + str(ndata_finite) +
                 ' ;' + str(ndata_S_N_gte_5) +
                 ' ;' + str(ndata_S_N_lt_0))
        n, bins, patches = plt.hist(xdata, bins=30,
                                    histtype='step',
                                    linewidth=linewidth,
                                    linestyle='--',
                                    range=xrange, label=label)


    plt.yscale('log')
    plt.xlim(xrange)

    plt.grid()
    plt.legend()

    plt.title(f'{fieldname}  \n All ;not NaN or Masked ;SN>=5 ;SN<0')
    plt.xlabel(fluxtype + ' _S_N')
    plt.ylabel('Number per bin')

    plotfile = 'DP1_' + fieldname + '_hist_' + fluxtype + '_S_N' + '.png'
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    logger.info('\n')
    plt.show()

    return


def plot_hist_redshift(table=None,
                       colname_redshift=None,
                       bins=30, zrange=(0.0, 6.0),
                       figsize=(8,6),
                       plot_title='',
                       plotfile_prefix=''):

    logger.info('\n')

    plt.figure(figsize=figsize)

    xcolname = colname_redshift

    xdata= table[xcolname]
    ndata = len(xdata)
    label = str(ndata)
    logger.info(f'{zrange} {bins}')
    plt.hist(xdata, bins=bins,
             histtype='step',
             color='black',
             linewidth=2,
             range=zrange,
             label=label)

    itest = (table['refExtendedness'] == 0)
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': refExtendness = 0'
    plt.hist(xdata, bins=bins,
                 histtype='step',
                 color='blue',
                 linewidth=2,
                 range=zrange,
                 label=label)

    itest = (table['refExtendedness'] != 0)
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': refExtendness != 0'
    plt.hist(xdata, bins=bins,
                 histtype='step',
                 color='red',
                 linewidth=2,
                 range=zrange,
                 label=label)

    plt.xlabel('Redshift')
    plt.ylabel('Number per bin')
    plt.title(plot_title)
    plt.legend()

    plotfile = plotfile_prefix + '_hist_redshift.png'
    logging.info(f'Saving plotfile: {plotfile}')
    plt.savefig(plotfile)

    #if showplots:
    logging.info('logging')
    plt.show()


    fig, axes = plt.subplots(1, 2, figsize=(12,6))

    itest = (table['refExtendedness'] == 0)
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': refExtendness = 0'
    axes[0].hist(xdata, bins=bins,
                 histtype='step',
                 color='blue',
                 linewidth=2,
                 range=zrange,
                 label=label)

    axes[0].legend()
    axes[0].set_xlabel('Redshift')
    axes[0].set_ylabel('Number per bin')



    itest = (table['refExtendedness'] != 0)
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': refExtendness != 0'
    axes[1].hist(xdata, bins=bins,
                 color='red',
                 linewidth=2,
                 histtype='step',
                 range=zrange,
                 label=label)

    axes[1].legend()
    axes[1].set_xlabel('Redshift')
    axes[1].set_ylabel('Number per bin')

    fig.suptitle(plot_title)

    plt.tight_layout()

    plotfile = plotfile_prefix + '_hist_redshift_byExtendedness.png'
    logging.info(f'Saving plotfile: {plotfile}')
    plt.savefig(plotfile)

    logger.info('\n')
    plt.show()

    return


def desi_plot_hist_redshift(table=None,
                            colname_redshift=None,
                            bins=30, zrange=(0.0, 6.0),
                            plot_title='',
                            plotfile_prefix=''):
    logger.info('\n')

    fig, axes = plt.subplots(1, 2, figsize=(12,6))

    xcolname = colname_redshift

    itest = (table['SPECTYPE'] == 'QSO')
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': QSO'
    axes[0].hist(xdata, bins=bins,
                 histtype='step',
                 color='black',
                 linewidth=2,
                 range=zrange,
                 label=label)

    itest = (table['SPECTYPE'] == 'QSO') & (table['refExtendedness'] == 0)
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': refExtendness = 0'
    axes[0].hist(xdata, bins=bins,
                 histtype='step',
                 color='blue',
                 linewidth=2,
                 range=zrange,
                 label=label)


    itest = (table['SPECTYPE'] == 'QSO') & (table['refExtendedness'] != 0)
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': refExtendness != 0'
    axes[0].hist(xdata, bins=bins,
                 histtype='step',
                 color='red',
                 linewidth=2,
                 range=zrange,
                 label=label)

    axes[0].legend()
    axes[0].set_xlabel('Redshift')
    axes[0].set_ylabel('Number per bin')


    itest = (table['SPECTYPE'] != 'QSO')
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': not QSO'
    axes[1].hist(xdata, bins=bins,
                 color='black',
                 linewidth=2,
                 histtype='step',
                 range=zrange,
                 label=label)


    itest = (table['SPECTYPE'] != 'QSO') & (table['refExtendedness'] == 0)
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': refExtendness = 0'
    axes[1].hist(xdata, bins=bins,
                 color='blue',
                 linewidth=2,
                 histtype='step',
                 range=zrange,
                 label=label)


    itest = (table['SPECTYPE'] != 'QSO') & (table['refExtendedness'] != 0)
    xdata= table[xcolname][itest]
    ndata = len(xdata)
    label = str(ndata) + ': refExtendness != 0'
    axes[1].hist(xdata, bins=bins,
                 color='red',
                 linewidth=2,
                 histtype='step',
                 range=zrange,
                 label=label)

    axes[1].legend()
    axes[1].set_xlabel('Redshift')
    axes[1].set_ylabel('Number per bin')

    fig.suptitle(plot_title)

    plt.tight_layout()

    plotfile = plotfile_prefix + '_hist_redshift_bySpecType.png'
    logging.info(f'Saving plotfile: {plotfile}')
    plt.savefig(plotfile)

    logger.info('\n')
    plt.show()

    return


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


def rd_table(infile=None, infostats=True):
    """
    read in a table and add metadata

    """
    print(f'Reading: {infile}')
    table = Table.read(infile)
    table.meta['Filename'] = infile
    logger.info(f'Table read in: {infile}')

    logger.info(f"Table: {table.meta['Filename']}")

    if infostats:
        table.info(['attributes', 'stats'])
        # logger.info(f'Elapsed time(secs): {time.time() - t0}\n')
        logger.info('\n')

    logger.info(f'Number of rows: {len(table)}')
    logger.info(f'Number of columns: ' +
                 f'{len(table.colnames)} {len(table.columns)}')
    # logger.info(f'Elapsed time(secs):  {time.time() - t0}\n')

    return table


import subprocess

def print_git_hash(short_hash=True):
    """
    Read and print the current git commit hash.

    Returns:
        str: The git commit hash, or None if not in a git repository

    If you want just the short hash (first 7 characters), you can change the git command to ['git', 'rev-parse', '--short', 'HEAD']
    """

    if short_hash:
        command = ['git', 'rev-parse', '--short', 'HEAD']

    if not short_hash:
        command = ['git', 'rev-parse', 'HEAD']
    try:
        # Run git command to get the current commit hash
        result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True)

        git_hash = result.stdout.strip()
        print(f"Git Hash: {git_hash}")
        return git_hash

    except (subprocess.CalledProcessError, FileNotFoundError):
        # Silently return None if not in a git repo or git not installed
        return None


# do some tests here
if __name__ == "__main__":

    githash_value = print_git_hash()
    if githash_value is None:
        print("No git hash available")

    githash_value = print_git_hash(short_hash=False)
    if githash_value is None:
        print("No git hash available")


    # define debug
    #debug = locals().get('debug', False)
    #print('debug: ', debug)

    args = parse_arguments()

    # define debug and DEBUG
    DEBUG = False
    debug = DEBUG
    if args.debug:
        DEBUG = True
        print(f"Debug mode is ON")

    # both for backward compatibility

    showplot = True
    showplots = showplot

    fieldnames = get_dp1_fieldnames()
    print(f'fieldnames: {fieldnames}')


    test_TempleModels = locals().get('debug', False)
    if test_TempleModels:
        # read the Temple models
        import TempleModels
        if debug:
            help(TempleModels)



        templemodels_dict = tm.get_colors_nested_dict()

        print()
        print(type(templemodels_dict))
        print(type(templemodels_dict['colors_M20_dict']))

        colors_list= list(templemodels_dict.keys())
        print(colors_list)

        for ikey, key in enumerate(templemodels_dict.keys()):
            print(ikey, key)


        print()
        print(f'showplots: {showplots}')
        print()


        colorcolorList = (['uggr', 'g', 'r', 'u', 'g', -1.5, 5.0, -1.5, 5.0],
                            ['uggr', 'u', 'g', 'g', 'r', -1.5, 5.0, -1.5, 5.0],
                            ['grri', 'r', 'i', 'g', 'r', -1.5, 5.0, -1.5, 5.0],
                            ['grri', 'g', 'r', 'r', 'i', -1.5, 5.0, -1.5, 5.0],
                            ['grrz', 'r', 'z', 'g', 'r', -1.5, 5.0, -1.5, 5.0],
                            ['grrz', 'g', 'r', 'r', 'z', -1.5, 5.0, -1.5, 5.0],
                            ['giiz', 'i', 'z', 'g', 'i', -1.5, 5.0, -1.5, 5.0],
                            ['giiz', 'g', 'i', 'i', 'z', -1.5, 5.0, -1.5, 5.0],
                            ['riiz', 'i', 'z', 'r', 'i', -1.5, 5.0, -1.5, 5.0],
                            ['izzY', 'z', 'Y', 'i', 'z', -1.5, 5.0, -1.5, 5.0],
                            ['zYYJ', 'Y', 'J', 'z', 'Y', -1.5, 5.0, -1.5, 5.0],
                            ['rzzW1', 'z', 'W1', 'r', 'z', -1.5, 5.0, -1.5, 5.0],
                            ['gzzW1', 'z', 'W1', 'g', 'z', -1.5, 5.0, -1.5, 5.0],
                            ['giW1W2', 'W1', 'W2', 'g', 'i',
                             -1.0, 2.0, -1.5, 5.0],
                            ['zW1W1W2', 'W1', 'W2', 'z', 'W1',
                             -1.0, 2.0, -1.5, 5.0])


        plot_redshiftrange = [0.0, 7.0]
        print(f'plot_redshiftrange: {plot_redshiftrange}')

        colormagsx = ['u', 'g']
        colormagsy = ['g', 'r']
        plot_colorrangex = (-1.5, 5.0)
        plot_colorrangey = (-1.5, 5.0)

        """
        tm.plot_color_color(templemodels_dict=templemodels_dict,
                           colormagsy=colormagsy,
                           absmaglist=
                           ['M20', 'M22', 'M24', 'M26', 'M28'],
                           plot_colorrangex=plot_colorrangex,
                           plot_colorrangey=plot_colorrangey)
        #                   showplots=showplots)
        """
