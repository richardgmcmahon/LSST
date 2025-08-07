import os
import sys
import time
import logging
import configparser

import numpy as np
import matplotlib.pyplot as plt

import TempleModels

from TempleModels import *
# help(TempleModels)

# help(TempleModels.plot_ugr)

from lsst.utils.plotting import (get_multiband_plot_colors,
                                 get_multiband_plot_symbols,
                                 get_multiband_plot_linestyles)

filter_colors = get_multiband_plot_colors()
filter_symbols = get_multiband_plot_symbols()
filter_linestyles = get_multiband_plot_linestyles()

field_filters = ['u', 'g', 'r', 'i', 'z', 'y']
bands = field_filters


logger = logging.getLogger(__name__)

def get_dp1_fieldnames():

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

     index_lsst = 1
     index_MQ = 1

     return index_lsst, index_MQ



def count_refband(table=None):
    """

    """
    unique, unique_counts = \
        np.unique(table['refBand'], return_counts=True)

    print(unique)
    print(unique_counts)
    for index, row in enumerate(unique):
        print(index, row, unique_counts[index])


    return

def plot_redshift_mag(table=None,
                      magtype='psfMag',
                      zrange=(0.0, 6.0),
                      markercolor=None,
                      markersize=4.0,
                      colname_redshift='z',
                      plot_title=None,
                      plotfile_prefix='',
                      savefig=True):
    """

    """
    # plot redshift v mag hubble diagram

    bands = ['u', 'g','r', 'i', 'z', 'y']
    magtype = 'psfMag'

    if plot_title is None:
        plot_title = table.meta['Filename']

    xcolname = colname_redshift
    for band in bands[0:6]:
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
        print(table[band + '_extendedness'])
        itest = (table[band + '_extendedness'] == 0)
        xdata = xdata_save[itest]
        ydata = ydata_save[itest]
        label = str(len(xdata)) + ' Pointlike'
        plt.plot(xdata, ydata, '.b', label=label)

        plt.xlim(zrange)
        plt.ylim(16.0, 28.0)

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


def plot_cmodel_psf(table=None,
                    band=None,
                    band_flaglist=None,
                    refExtendedness=False,
                    markersize=4.0,
                    xrange=None,
                    yrange=None,
                    multiplot=False,
                    plot_title=None,
                    plotfile_prefix='',
                    savefig=True):
    """

    plot cmodel - psf


    """

    # global logger

    logging.info('')

    if plot_title is None:
        plot_title = table.meta['Filename']

    if band is None:
        bands = ['u', 'g','r', 'i', 'z', 'y']
    if band is not None:
        bands = band

    magtype1 = 'psfMag'
    magtype2 = 'cModelMag'



    table['refBand'].info(['attributes', 'stats'])

    unique, unique_counts = \
            np.unique(table['refBand'], return_counts=True)
    print('refBand', len(table))
    print(unique)
    print(unique_counts)
    for index, row in enumerate(unique):
        print(index, row, unique_counts[index])

    ndata_finite = np.count_nonzero(
            ~np.isnan(table['refExtendedness']))
    table['refExtendedness'].info(['attributes', 'stats'])
    print(f'Number of finite refExtendedness values : {ndata_finite}')

    table_save = table
    for band in bands[0:6]:

        table = table_save

        if band_flaglist is not None:
            for band_flag in band_flaglist:
                itest = (table[band + '_' + band_flag] == 0)
                table = table[itest]

        print()
        logging.info('')
        plot_suptitle = band + '_extendedness'
        ndata = len(table)
        print(table[band + '_extendedness'])

        ndata_finite = np.count_nonzero(
            ~np.isnan(table[band + '_extendedness']))
        table[band + '_extendedness'].info(['attributes', 'stats'])
        print(f'Number of finite {band} ExtendedNess values : {ndata_finite}')
        print()

        ndata_finite = np.count_nonzero(
            ~np.isnan(table[band + '_psfMag']))
        table[band + '_psfMag'].info(['attributes', 'stats'])
        print(f'Number of finite {band}_psfMag values : {ndata_finite}')
        print()

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

        unique, unique_counts = \
            np.unique(table['refExtendedness'], return_counts=True)
        print('refExtendedness')
        print(unique)
        print(unique_counts)
        for index, row in enumerate(unique):
            print(index, row, unique_counts[index])
        print()

        label = str(ndata_finite) + '/' + str(ndata)
        plt.plot(np.nan, np.nan, color='None', label=label)

        xlabel = band + '_' + magtype1 + ' - ' + \
            band + '_' + magtype2
        xdata_save = table[band + '_' + magtype1] - \
            table[band + '_' + magtype2]

        ndata_finite = np.count_nonzero(~np.isnan(xdata_save))
        print(f'Number of finite {band}_{magtype1} - {band}_{magtype2} values : {ndata_finite}')


        ycolname = band + '_' + magtype1
        ydata_save = table[ycolname]

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
            plot_suptitle = 'refExtendedness'

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

        plt.legend(loc='upper right')
        plt.grid()
        plt.xlabel(xlabel)
        plt.ylabel(ycolname)
        plt.title(plot_title)
        plt.gca().invert_yaxis()
        plt.suptitle(plot_suptitle)

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
                         plot_title=None,
                         plotfile_prefix='',
                         savefig=True):
    # redshift versus color

    if plot_title is None:
        plot_title = table.meta['Filename']

    xcolname = colname_redshift

    for iband, band in enumerate(bands[0:5]):

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

        print(table[band + '_extendedness'])
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
        plt.ylabel(ycolname)

        logging.info('')
        plt.show()

    return


def plot_color_mag(
        table=None,
        survey=None,
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

    filename = table.meta['Filename']
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
    footnote = timestamp + ': ' + \
        os.path.basename(__file__) + ': ' + \
        plotfile
    plt.figtext(0.01, 0.01, footnote,
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


    filename = table.meta['Filename']
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
    footnote = timestamp + ': ' + \
        os.path.basename(__file__) + ': ' + \
        plotfile
    plt.figtext(0.01, 0.01, footnote,
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

    help(plt.hexbin)

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

    plt.show()

    return


def plot_image_quality():
    """


    """

    return



def stats_column(table=None, colname=None):


    logger.info('\n')
    table[colname].info(['attributes', 'stats'])
    data = table[colname]

    is_masked = np.ma.is_masked(data)
    print(f'is_masked: {is_masked}')

    if is_masked:
        masked_data = data[data.mask]
        unmasked_data = data[~data.mask]
        print('Number of rows; all, unmasked, masked: ' +
          f'{len(data)} {len(unmasked_data)} {len(data_masked)}')

    print(type(data))
    print(f'data[0, -1]: {data[0]} {data[-1]}')

    ndata_finite = np.count_nonzero(np.isfinite(data))
    print(f'Number of finite values: {ndata_finite}')

    print()
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
             '/' + str(ndata_finite) +
             '/' + str(ndata_S_N_gte_5) +
             '/' + str(ndata_S_N_lt_0))
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
                 '/' + str(ndata_finite) +
                 '/' + str(ndata_S_N_gte_5) +
                 '/' + str(ndata_S_N_lt_0))
        n, bins, patches = plt.hist(xdata, bins=30,
                                    histtype='step',
                                    linewidth=linewidth,
                                    linestyle='--',
                                    range=xrange, label=label)


    plt.yscale('log')
    plt.xlim(xrange)

    plt.grid()
    plt.legend()

    plt.title(f'{fieldname}  \n All \\ not NaN or Masked\\ SN>=5\\ SN<0')
    plt.xlabel(fluxtype + ' _S_N')
    plt.ylabel('Number per bin')

    plotfile = 'DP1_' + fieldname + '_' + fluxtype + '_S_N' + '.png'
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    logging.info('logging')
    logger.info('logger')

    plt.show()

    return



# do some tests here
if __name__ == "__main__":

    # read the Temple models
    import TempleModels
    help(TempleModels)

    import os
    import sys
    import time


    fieldnames = get_dp1_fieldnames()
    print(f'fieldnames: {fieldnames}')


    showplots = True

    templemodels_dict = get_colors_nested_dict()

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
    plot_color_color(templemodels_dict=templemodels_dict,
                       colormagsy=colormagsy,
                       absmaglist=
                       ['M20', 'M22', 'M24', 'M26', 'M28'],
                       plot_colorrangex=plot_colorrangex,
                       plot_colorrangey=plot_colorrangey)
    #                   showplots=showplots)
