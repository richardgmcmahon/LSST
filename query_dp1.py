"""

https://cds.unistra.fr/help/documentation/vizier-more/adql-vizier/

https://sdm-schemas.lsst.io/dp1.html

https://sdm-schemas.lsst.io/dp1.html#Object

https://sdm-schemas.lsst.io/v/

https://rtn-095.lsst.io/

https://rtn-095.lsst.io/v/

https://dp1.lsst.io/

https://pyvo.readthedocs.io/

https://dp1.lsst.io/tutorials/notebook/301/notebook-301-0.html

https://dp1.lsst.io/tutorials/notebook/301/notebook-301-4.html

https://github.com/rubin-dp0/tutorial-notebooks/blob/main/DP02_13a_Image_Cutout_SciDemo.ipynb

SELECT COUNT(*) FROM dp1.CoaddPatches

SELECT TOP 100 * FROM dp1.CoaddPatches

Visits
855

=> LSSTCam is 21 rafts, so this is like 855/21 = 41  LSSTCam visits



coaddPatches
652

SELECT COUNT(DISTINCT lsst_tracts) FROM FROM dp1.CoaddPatches
29

SELECT lsst_tract, COUNT(lsst_tract) AS Number
FROM dp1.CoaddPatches
GROUP BY lsst_tract


Object
SELECT COUNT(*) FROM FROM dp1.Object
2299757

Source
45,565,632

ForcedSource
268,796,943



https://community.lsst.org/t/how-to-obtain-local-access-to-dp0-2-images/8540

"""

import pyvo
import os
import sys
import time

#
import argparse
import configparser
import logging
import getpass
import inspect
import traceback


from astropy.table import Table
import astropy.units as u
import numpy as np

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

#from lsst.rsp import get_tap_service

# needs pip install lsst_utils
# conda install wrapt
from lsst.utils.plotting import (get_multiband_plot_colors,
                                 get_multiband_plot_symbols,
                                 get_multiband_plot_linestyles)


# my functions
import lsst_util as lu

filter_colors = get_multiband_plot_colors()
filter_symbols = get_multiband_plot_symbols()
filter_linestyles = get_multiband_plot_linestyles()


field_filters = ['u', 'g', 'r', 'i', 'z', 'y']
bands = field_filters


def connect_rsp_tap():
    """
    https://dp0-2.lsst.io/data-access-analysis-tools/api-intro.html#use-of-pyvo-with-the-rsp-tap-service

    https://community.lsst.org/t/url-to-access-images-via-sia/10480

    """

    logger.info('')

    RSP_TAP_SERVICE = 'https://data.lsst.cloud/api/tap'
    homedir = os.path.expanduser('~')
    token_file = os.path.join(homedir,'.rsp-tap.token')
    with open(token_file, 'r') as f:
        token_str = f.readline()

    print(f'token: {token_str}')
    cred = pyvo.auth.CredentialStore()
    cred.set_password("x-oauth-basic", token_str)
    credential = cred.get("ivo://ivoa.net/sso#BasicAA")
    rsp_tap = pyvo.dal.TAPService(RSP_TAP_SERVICE, session=credential)

    return rsp_tap


def run_async_query(service=None, query=None):
    """


    """
    rsp_tap = service

    print('Run async tap query:', )
    logger.info('')
    #logger.info(f'{job.url}')

    job = rsp_tap.submit_job(query)
    # help(job.run)
    job.run()
    job.wait(phases=['COMPLETED', 'ERROR'])
    print('Job phase is', job.phase)
    if job.phase == 'ERROR':
        job.raise_if_error()

    logger.info('')
    print('Elapsed time(secs): ',time.time() - t0)
    assert job.phase == 'COMPLETED'

    logger.info('Fetch results table')
    results = job.fetch_result().to_table()
    print(f"Total number of rows: {len(results)}")

    logger.info('')
    # results= rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)

    # objtab = results.to_table()

    print()
    results.info(['attributes', 'stats'])
    logger.info('')

    return results


def query_schema(rsp_tap=None):
    """

    SELECT * FROM tap_schema.schemas

    SELECT schema_name, table_name, description
    FROM tap_schema.tables
    WHERE schema_name != 'tap_schema'
    ORDER BY schema_name, table_name


    """

    query_number = 2

    query1 = "SELECT * FROM tap_schema.schemas"

    query2 = """
    SELECT schema_name, table_name, description
    FROM tap_schema.tables
    WHERE schema_name != 'tap_schema'
    ORDER BY schema_name, table_name """

    if query_number == 1:
        query = query1

    if query_number == 2:
        query = query2

    results = rsp_tap.run_sync(query)

    table = results.to_table()

    table.info(['attributes', 'stats'])

    if query_number == 1:
        print(table['schema_name'])

    if query_number == 2:
        print(table['table_name'])

    return table


def query_CoaddPatches(top=None):
    """CoaddPatches


    lsst_tract (ID number of the top level, 'tract',
    within the standard LSST skymap)
    lsst_patch (ID number of the second level, 'patch',
    within the standard LSST skymap)
    s_ra (Central Spatial Position in ICRS; Right ascension)
    s_dec (Central Spatial Position in ICRS; Declination)
    s_region (Sky region covered by the coadd (expressed in ICRS frame))
    e.g. POLYGON ICRS 6.270469 -72.819733 6.909381 -72.813516 6.885285 -72.624796 6.253075 -72.630941

    """
    # Count the rows
    query = "SELECT COUNT(*) FROM dp1.CoaddPatches"
    print(query)
    results = rsp_tap.run_sync(query)
    results = results.to_table()
    print(f'Number of rows: {results[0][0]}')
    print()


    # get the column names
    query = "SELECT TOP 10 * FROM dp1.CoaddPatches"
    print(query)
    results = rsp_tap.run_sync(query)
    results = results.to_table()
    results.info(['attributes'])


    # Distinct tracts
    query_lines = ["SELECT "]
    if top is not None:
        query_lines.append(f"       TOP {top} ")
    query_lines.append(f"       DISTINCT lsst_tract ")
    query_lines.append(f"       FROM dp1.CoaddPatches")
    if all == False:
        query_lines.append(
            f"WHERE CONTAINS(POINT('ICRS', s_ra, s_dec), "
            f"CIRCLE('ICRS', {ra_centre}, {dec_centre}, {radius})) = 1 "
        )
    query_lines.append("        ORDER BY lsst_tract ASC")

    query = "\n".join(query_lines)

    print(query)
    print()
    results = rsp_tap.run_sync(query)
    results = results.to_table()
    print(results)
    results.info(['attributes', 'stats'])


    # Distinct patchs
    query_lines = ["SELECT "]
    if top is not None:
        query_lines.append(f"       TOP {top} ")
    query_lines.append(f"       DISTINCT lsst_tract, lsst_patch")
    query_lines.append(f"       FROM dp1.CoaddPatches")
    if all == False:
        query_lines.append(
            f"WHERE CONTAINS(POINT('ICRS', s_ra, s_dec), "
            f"CIRCLE('ICRS', {ra_centre}, {dec_centre}, {radius})) = 1 "
        )
    query_lines.append("        ORDER BY lsst_tract ASC")

    query = "\n".join(query_lines)

    print(query)
    print()

    results = rsp_tap.run_sync(query)

    results = results.to_table()

    print(results)
    results.info(['attributes', 'stats'])


    query = """
            SELECT DISTINCT lsst_tract
            FROM dp1.CoaddPatches
            WHERE CONTAINS(POINT('ICRS', s_ra, s_dec), CIRCLE('ICRS', {}, {}, {}))=1
            ORDER BY lsst_tract
            """.format(ra_centre, dec_centre, radius)

    print(f'query: {query}')

    results = rsp_tap.run_sync(query)

    results = results.to_table()

    print(results)
    results.info(['attributes', 'stats'])

    tracts = results['lsst_tract']



    query = """
            SELECT DISTINCT lsst_tract, lsst_patch
            FROM dp1.CoaddPatches
            WHERE CONTAINS(POINT('ICRS', s_ra, s_dec), CIRCLE('ICRS', {}, {}, {}))=1
            ORDER BY lsst_tract
            """.format(ra_centre, dec_centre, radius)

    print(f'query: {query}')

    results = rsp_tap.run_sync(query)

    results = results.to_table()

    print(results)
    results.info(['attributes', 'stats'])

    tracts = results['lsst_tract']


    return tracts



def query_visits():
    """


    """


    query = """
            SELECT visit, band, expMidptMJD, expMidpt
            FROM dp1.Visit
            WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {}, {}, {}))=1
            ORDER BY expMidptMJD
            """.format(ra_centre, dec_centre, radius)

    print()
    print(query)
    print()
    t0 = time.time()
    results = rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)

    visits = results.to_table()

    print()
    visits.info(['attributes', 'stats'])
    print(f"Total number of visits: {len(visits)}")

    all_bands = np.array([visit['band'] for visit in visits], dtype=str)
    unique_bands, counts = np.unique(all_bands, return_counts=True)

    for band, count in zip(unique_bands, counts):
        print(band, count)

    return visits


def query_ccdvisits(plotfile_prefix=None):
    """CcdVisits
    """



    title = ''
    if plotfile_prefix is not None:
        title = plotfile_prefix
        plotfile_prefix = plotfile_prefix + '_'


    query = """
            SELECT visitId, ra, dec, band, magLim, seeing
            FROM dp1.CcdVisit
            WHERE CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {}, {}, {}))=1
            ORDER BY visitId
            """.format(ra_cen, dec_cen, radius)

    results = rsp_tap.run_sync(query)

    ccd_visits = results.to_table()

    print()
    ccd_visits.info(['attributes', 'stats'])

    for band in field_filters:
        print(f"Median in {band}-band: {np.ma.median(ccd_visits[ccd_visits['band'] == band]['seeing']): .2f}")


    min_seeing = np.min(ccd_visits['seeing'])
    max_seeing = np.max(ccd_visits['seeing'])
    seeing_bins = np.arange(min_seeing, max_seeing, 0.05)

    fig = plt.figure(figsize=(7, 5))
    for band in field_filters:
        print(f"Median in {band}-band: \
                {np.ma.median(ccd_visits[ccd_visits['band'] == band]['seeing']): .2f}")
        label = band + '/' + str(len(ccd_visits[ccd_visits['band'] == band]))
        n, bins, patches = plt.hist(ccd_visits[ccd_visits['band'] == band]['seeing'],
                                    seeing_bins, histtype='step',
                                    linewidth=2,
                                    color=filter_colors[band],
                                    label=label)
        for patch in patches:
            patch.set_linestyle(filter_linestyles[band])

    plt.legend(loc='upper right')
    plt.xlabel('PSF FWHM (arcsec)')
    plt.ylabel('Number of visits')
    plt.minorticks_on()
    plotfile = plotfile_prefix + 'ccdvisits_psf_fwhm' + '.png'
    print('Saving:', plotfile)
    plt.savefig(plotfile)
    plt.show()


    min_mag = np.min(ccd_visits['magLim'])
    max_mag = np.max(ccd_visits['magLim'])
    maglim_bins = np.arange(min_mag, max_mag, 0.1)

    fig = plt.figure(figsize=(7, 5))

    for band in field_filters:
        label = band + '/' + \
            str(len(ccd_visits[ccd_visits['band'] == band]['magLim'])) + '/' + \
            str(len(ccd_visits['band']))

        print(f"Median in {band}-band: \
                {np.ma.median(ccd_visits[ccd_visits['band'] == band]['magLim']): .2f}")
        n, bins, patches = plt.hist(ccd_visits[ccd_visits['band'] ==
                                               band]['magLim'],
                                    maglim_bins, histtype='step',
                                    linewidth=2,
                                    color=filter_colors[band],
                                    label=label)
        for patch in patches:
            patch.set_linestyle(filter_linestyles[band])

    plt.legend(loc='upper left')
    plt.title(title)
    plt.xlabel('limiting magnitude')
    plt.ylabel('Number of visits')
    plt.minorticks_on()

    plotfile = plotfile_prefix + 'ccdvisits_nvisits' + '.png'
    print('Saving:', plotfile)
    plt.savefig(plotfile)
    plt.show()


    plt.show()

    return


def query_Object(rsp_tab=None):
    """


    4.1. Magnitude distributions

    Plot histograms of object magnitudes, separating the samples into
    stars and galaxies.

    Use the refExtendedness flag to distinguish them: treat objects with
    refExtendedness == 1 as likely galaxies (extended) and those with
    refExtendedness == 0 as likely stars (point sources).

    Use cModelMag mags for galaxies and psfMag for stars.

    SELECT tract, COUNT(*) AS Number
    FROM dp1.Object
    GROUP BY tract

    SELECT tract, patch COUNT(*) AS Number
    FROM dp1.Object
    GROUP BY tract, patch
    608


    """

    query = """SELECT tract, patch, COUNT(*) AS Number
           FROM dp1.Object
           GROUP BY tract, patch """

    print(f'Query: {query}')
    t0 = time.time()
    results= rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)
    results = results.to_table()
    results.info(['attributes', 'stats'])
    print(results)

    Number = results['Number']
    nPatches = len(results)
    print(f'Number of Patches: {nPatches}')
    print(f'Number of Sources per Patch [min, max, median, 10%, 90%]: ' +
          f'{np.min(Number)}, ' +
          f'{np.max(Number)}, {int(np.median(Number))}, ' +
          f'{int(np.percentile(Number, 10))}, ' +
          f'{int(np.percentile(Number, 90))}')


    query = """SELECT tract, COUNT(*) AS Number
           FROM dp1.Object
           GROUP BY tract"""

    print(f'Query: {query}')
    t0 = time.time()
    results= rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)
    results = results.to_table()
    results.info(['attributes', 'stats'])
    print(results)

    Number = results['Number']
    nTracts = len(results)
    print(f'Number of Tracts: {nTracts}')
    print(f'Number of Sources per Tract [min, max, median, 10%, 90%]: ' +
          f'{np.min(Number)}, ' +
          f'{np.max(Number)}, {int(np.median(Number))}, ' +
          f'{int(np.percentile(Number, 10))}, ' +
          f'{int(np.percentile(Number, 90))}')



    return

    query = """SELECT objectId, coord_ra, coord_dec, ebv,
            {}_psfMag, {}_cModelMag, {}_psfMag, {}_cModelMag,
            {}_psfMag, {}_cModelMag, {}_psfMag, {}_cModelMag,
            {}_psfMag, {}_cModelMag, {}_psfMag, {}_cModelMag,
            refExtendedness
            FROM dp1.Object
            WHERE CONTAINS(POINT('ICRS', coord_ra, coord_dec),
                  CIRCLE('ICRS', {}, {}, {})) = 1
            ORDER BY objectId ASC
            """.format(field_filters[0], field_filters[0],
                       field_filters[1], field_filters[1],
                       field_filters[2], field_filters[2],
                       field_filters[3], field_filters[3],
                       field_filters[4], field_filters[4],
                       field_filters[5], field_filters[5],
                       ra_centre, dec_centre, radius)
    print(query)

    t0 = time.time()
    results= rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)

    objtab = results.to_table()

    print()
    objtab.info(['attributes', 'stats'])

    outfile = 'objtab_dp1_' + fieldname + '.fits'
    objtab.write(outfile, overwrite=True)

    ptsource = (objtab['refExtendedness'] == 0)

    mag_bins = np.arange(15.4, 28.6, 0.2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    plt.subplots_adjust(hspace=0, wspace=0)

    for band in field_filters:
        mag_psf = objtab[f"{band}_psfMag"]
        mag_cmodel = objtab[f"{band}_cModelMag"]

        valid_mags = (15 < mag_psf) & (mag_psf < 30) & (15 < mag_cmodel) & (mag_cmodel < 30)

        print(f"Number of objects in {band}-band: {np.sum(valid_mags)}")

        n_psf, bins_psf, patches_psf = ax1.hist(mag_psf[valid_mags & ptsource], bins=mag_bins,
                                                histtype='step', linewidth=2,
                                                color=filter_colors[band], label=band,)
        for patch in patches_psf:
            patch.set_linestyle(filter_linestyles[band])

        n_cmodel, bins_cmodel, patches_cmodel = ax2.hist(
            mag_cmodel[valid_mags & ~ptsource],
            bins=mag_bins,
            histtype='step',
            linewidth=2,
            color=filter_colors[band],
            label=band,
        )
        for patch in patches_cmodel:
            patch.set_linestyle(filter_linestyles[band])

    ax1.set_xlim(mag_bins.min(), mag_bins.max())
    ax1.set_yscale('log')
    ax1.set_xlabel('Magnitude')
    ax1.set_ylabel('Number of objects')
    ax1.set_title('Point sources (PSF mags)')
    ax1.legend(loc='upper left', ncols=2)
    ax1.minorticks_on()

    ax2.set_xlim(mag_bins.min(), mag_bins.max())
    ax2.set_xlabel('Magnitude')
    ax2.set_title('Extended sources (cModel mags)')
    ax2.minorticks_on()

    plotfile = plotfile_prefix + 'Object_hist_mag' + '.png'
    print('Saving:', plotfile)
    plt.savefig(plotfile)
    plt.show()


    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        2, 2, figsize=(7, 6), height_ratios=[1.5, 2.5])
    plt.subplots_adjust(hspace=0, wspace=0)

    grmin, grmax = -0.9, 2.3
    rimin, rimax = -1.3, 2.8
    magmin, magmax = 15.8, 26.8

    ax1.hexbin(objtab[ptsource]['g_psfMag']-objtab[ptsource]['r_psfMag'],
               objtab[ptsource]['r_psfMag']-objtab[ptsource]['i_psfMag'],
               gridsize=(150, 150),
               extent=(grmin, grmax, rimin, rimax), bins='log', cmap='Grays')
    ax1.set_title('Point sources (PSF mags)')
    ax1.set_xlim(grmin, grmax)
    ax1.set_ylim(rimin, rimax)
    ax1.set_ylabel(r'$(r-i)$')
    ax1.set_xticklabels([])
    ax1.minorticks_on()

    ax2.hexbin(objtab[~ptsource]['g_cModelMag']-objtab[~ptsource]['r_cModelMag'],
               objtab[~ptsource]['r_cModelMag']-objtab[~ptsource]['i_cModelMag'],
               gridsize=(150, 150),
               extent=(grmin, grmax, rimin, rimax), bins='log', cmap='Grays')
    ax2.set_title('Extended sources (cModel mags)')
    ax2.set_xlim(grmin, grmax)
    ax2.set_ylim(rimin, rimax)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.minorticks_on()

    ax3.hexbin(objtab[ptsource]['g_psfMag']-objtab[ptsource]['r_psfMag'],
               objtab[ptsource]['r_psfMag'], gridsize=(150, 200),
               extent=(grmin, grmax, magmin, magmax), bins='log', cmap='Grays')

    ax3.set_xlim(grmin, grmax)
    ax3.set_ylim(magmax, magmin)
    ax3.set_xlabel(r'$(g-r)$')
    ax3.set_ylabel(r'$r$ magnitude')
    ax3.minorticks_on()

    ax4.hexbin(objtab[~ptsource]['g_cModelMag']-objtab[~ptsource]['r_cModelMag'],
               objtab[~ptsource]['r_cModelMag'], gridsize=(150, 200),
               extent=(grmin, grmax, magmin, magmax), bins='log', cmap='Grays')
    ax4.set_xlim(grmin, grmax)
    ax4.set_ylim(magmax, magmin)
    ax4.set_xlabel(r'$(g-r)$')
    ax4.set_yticklabels([])
    ax4.minorticks_on()

    plt.show()

    return objtab


def debug_datalink_url(url):
    import requests

    homedir = os.path.expanduser('~')
    token_file = os.path.join(homedir,'.rsp-tap.token')
    with open(token_file, 'r') as f:
        token_str = f.readline()

    cred = pyvo.auth.CredentialStore()
    cred.set_password("x-oauth-basic", token_str)
    credential = cred.get("ivo://ivoa.net/sso#BasicAA")
    print(f'credential: {credential}')

    # Fetch the content at the URL
    print(f"Fetching: {url}")
    response = requests.get(url, auth=credential)
    print(f"HTTP status code: {response.status_code}")
    print(f"Content-Type: {response.headers.get('Content-Type')}")

    # Save the content locally for inspection
    with open('debug_votable.xml', 'wb') as f:
        f.write(response.content)
    print("Saved response to 'debug_votable.xml'")

    # Print first few lines to inspect if it looks like a VOTABLE
    print("First 500 characters of response:")
    print(response.text[:500])

    # Now try to parse with PyVO, but catch exceptions
    try:
        dl_results = DatalinkResults.from_result_url(url)
        print("Parsed VOTABLE successfully!")
        return dl_results
    except Exception as e:
        print("Exception occurred while parsing VOTABLE:")
        print(e)



def query_ObsCore():
    """
    https://dp1.lsst.io/tutorials/notebook/103/notebook-103-3.html

    ivoa.ObsCore
    dataproduct_type (Data product (file content) primary type)
    dataproduct_subtype (Data product specific type)
    calib_level (Calibration level of the observation: in {0, 1, 2, 3, 4})
    lsst_band (Abstract filter band designation)
    em_min (start in spectral coordinates)
    em_max (stop in spectral coordinates)
    lsst_tract (Upper level of LSST coadd skymap hierarchy)
    lsst_patch (Lower level of LSST coadd skymap hierarchy)
    lsst_filter (Physical filter designation from the LSSTCam filter set)
    lsst_visit (Identifier for a specific LSSTCam pointing)
    lsst_detector (Identifier for CCD within the LSSTCam focal plane)
    t_exptime (Total exposure time)
    t_min (Start time in MJD)
    t_max (Stop time in MJD)
    s_ra (Central Spatial Position in ICRS; Right ascension)
    s_dec (Central Spatial Position in ICRS; Declination)
    s_fov (Estimated size of the covered region as the diameter of a containing circle)
    obs_id (Internal ID given by the ObsTAP service)
    obs_collection (Name of the data collection)
    o_ucd (Nature of the observable axis)
    facility_name (The name of the facility, telescope, or space craft used for the observation)
    instrument_name (The name of the instrument used for the observation)
    obs_title (Brief description of dataset in free format)
    s_region (Sky region covered by the data product (expressed in ICRS frame))
    access_url (URL used to access dataset)
    access_format (Content format of the dataset)
    obs_publisher_did (ID for the Dataset given by the publisher)
    target_name (Object of interest)
    s_resolution (Spatial resolution of data as FWHM of PSF)
    s_xel1 (Number of elements along the first coordinate of the spatial axis)
    s_xel2 (Number of elements along the second coordinate of the spatial axis)
    t_resolution (Temporal resolution FWHM)
    t_xel (Number of elements along the time axis)
    pol_xel (Number of elements along the polarization axis)
    em_xel (Number of elements along the spectral axis)
    em_res_power (Value of the resolving power along the spectral axis (R))


    target_ra = 53.2
    target_dec = -28.1


    """

    from pyvo.dal.adhoc import DatalinkResults

    from lsst.rsp import get_tap_service
    from lsst.rsp.utils import get_pyvo_auth

    logger.info('')

    # help(get_pyvo_auth)

    service = get_tap_service("tap")

    RSP_TAP_SERVICE = 'https://data.lsst.cloud/api/tap'
    homedir = os.path.expanduser('~')
    token_file = os.path.join(homedir,'.rsp-all.token')
    print(f'Reading: {token_file}')
    with open(token_file, 'r') as f:
        token_str = f.readline()

    print(f'token: {token_str}')
    cred = pyvo.auth.CredentialStore()
    cred.set_password("x-oauth-basic", token_str)
    credential = cred.get("ivo://ivoa.net/sso#BasicAA")
    print(f'credential: {credential}')


    service = pyvo.dal.TAPService(RSP_TAP_SERVICE, session=credential)
    auth_session = service._session
    print(auth_session)

    query1 = """
            SELECT * FROM ivoa.ObsCore
            WHERE CONTAINS(POINT('ICRS', 53.2,-28.1), s_region) = 1
            """

    query2 = """
            SELECT lsst_visit, lsst_band,
                   lsst_tract, lsst_patch,
                   s_ra, s_dec,
                   s_region,
                   access_url
            FROM ivoa.ObsCore
            WHERE (calib_level = 3)
            AND (lsst_band = 'y')
            AND (dataproduct_subtype = 'lsst.deep_coadd')
            AND CONTAINS(POINT('ICRS', {}, {}), s_region) = 1
            """.format(ra_centre, dec_centre, radius)
#            """.format(ra_cen, dec_cen, radius)

    query3 = """SELECT * FROM ivoa.ObsCore
        WHERE CONTAINS(POINT('ICRS', {},{}), s_region) = 1
        AND dataproduct_subtype='lsst.deep_coadd'
        AND lsst_band = 'i'
        ORDER BY lsst_tract
        """.format(target_ra, target_dec)

    query = query3

    print('query:', query)
    t0 = time.time()
    results = service.run_sync(query)

    #results = service.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)

    results = results.to_table()
    print()
    results.info(['attributes', 'stats'])
    print(f"Total number of images: {len(results)}")

    values, counts = np.unique(results['dataproduct_subtype'],
                           return_counts=True)
    for value, count in zip(values, counts):
        print(value, count)

    help(DatalinkResults)

    help(DatalinkResults.from_result_url)

    auth_session = service._session
    print(f'auth_session: {auth_session}')

    print(f'credential: {credential}')

    datalink_url = results['access_url'][0]
    print(f'datalink_url {datalink_url}\n')

    # debug_datalink_url(datalink_url)

    """
    https://github.com/rubin-dp0/tutorial-notebooks/blob/main/DP02_13a_Image_Cutout_SciDemo.ipynb
    """

    dataLinkUrl = datalink_url
    # dataLinkUrl = results[0].getdataurl()

    import requests
    response = requests.get(datalink_url)
    print(response.text)

    """
    import requests

    url = "https://example.com/protected/resource"
    headers = {"Authorization": "Bearer YOUR_TOKEN"}
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
       print("Authentication successful!")
    else:
        print(f"Authentication failed: {response.status_code}")
        print(response.text)


    headers = {"Authorization": token_str}
    response = requests.get(datalink_url, headers=headers)

    if response.status_code == 200:
       print("Authentication successful!")
    else:
        print(f"Authentication failed: {response.status_code}")
        print(response.text)

    """

    # OR

    """
    import requests
    from pyvo.dal.adhoc import DatalinkResults

    session = requests.Session()
    session.auth = ('your_username', 'your_password')  # if using basic auth
    # OR for tokens:
    session.headers.update({'Authorization': 'Bearer YOUR_TOKEN'})

    dl_results = DatalinkResults.from_result_url(result_url, session=session)

    """

    logger.info('')

    session = requests.Session()
    session.auth = ('x-oauth-basic',
                    'gt-51CZiN4dz25_ByHLl6G3PA.4roTRC8psE8Xg86wNDqjaw')  # if using basic auth

    auth_session = service._session
    dl_results = DatalinkResults.from_result_url(
         dataLinkUrl, session=auth_session)

    print(f"Datalink status: {dl_results.status}. Datalink service url: {dataLinkUrl}")



    #dl_result = DatalinkResults.from_result_url(datalink_url, session=get_pyvo_auth())
    #dl_result = DatalinkResults.from_result_url(
    #    datalink_url, session=auth_session)
    dl_result = DatalinkResults.from_result_url(
        datalink_url, session=credential)

    logger.info('')

    #import pyvo
    import requests
    token = "gt-HahpPOxt0KHlbtjKewA7aw.rJFmhejAwSh-FjGavXWP9g"
    session = requests.Session()
    session.headers["Authorization"] = f"Bearer {token}"
    dl_result = DatalinkResults.from_result_url(
        datalink_url, session=session)


    logger.info('')
    image_url = dl_result.getrecord(0).get('access_url')

    sys.exit()

    print(obscore['access_url'][0])

    from pyvo.dal.adhoc import DatalinkResults
    # https://pyvo.readthedocs.io/en/latest/dal/index.html#datalink
    # In this example you know the URL from somewhere

    url = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/datalink?ID=ivo%3A%2F%2Fcadc.nrc.ca%2FHSTHLA%3Fhst_12477_28_acs_wfc_f606w_01%2Fhst_12477_28_acs_wfc_f606w_01_drz'

    url = obscore['access_url'][0]

    # https://dp1.lsst.io/tutorials/notebook/103/notebook-103-3.html
    print(f'datalink: {url} \n')

    from lsst.rsp import get_tap_service
    from lsst.rsp.utils import get_pyvo_auth



    #service = get_tap_service("tap")


    service = rsp_tap
    query = "SELECT column_name, datatype, description, unit " \
        "FROM tap_schema.columns " \
        "WHERE table_name = 'ivoa.ObsCore'"
    print(query)
    results = service.search(query).to_table()
    print(results)


    dl_result = DatalinkResults.from_result_url(url, session=get_pyvo_auth())
    image_url = dl_result.getrecord(0).get('access_url')

    print(f'image_url: {image_url}')

    sys.exit()

    import requests
    f = requests.get(url)
    print(f.text)

    sys.exit()

    import urllib.request
    f = urllib.request.urlopen(url)
    myfile = f.read()
    print(myfile)

    sys.exit()


    help(DatalinkResults.from_result_url)

    homedir = os.path.expanduser('~')
    token_file = os.path.join(homedir,'.rsp-tap.token')
    with open(token_file, 'r') as f:
        token_str = f.readline()

    cred = pyvo.auth.CredentialStore()
    cred.set_password("x-oauth-basic", token_str)
    credential = cred.get("ivo://ivoa.net/sso#BasicAA")

    dl_results = DatalinkResults.from_result_url(
        url, session=credential)

    print(dl_results)

    next(datalink.bysemantics("#this")).content_type


    return


def select_highredshift_galaxies_grz():

    title = 'dp1_ecdfs'

    query_lines = [
        "SELECT objectId, coord_ra, coord_dec, tract, patch, "
]
    for filt in field_filters:
        query_lines.append(f"       {filt}_cModelFlux, {filt}_cModelFluxErr, ")
        query_lines.append(f"       {filt}_cModelMag, {filt}_cModelMagErr, ")
        query_lines.append(f"       {filt}_extendedness, ")

    query_lines.append("       refExtendedness, refBand")

    query_lines.append("FROM dp1.Object")

    query_lines.append(
        f"WHERE CONTAINS(POINT('ICRS', coord_ra, coord_dec), "
        f"CIRCLE('ICRS', {ra_cen}, {dec_cen}, {radius})) = 1 "
        #f"AND r_cModelFlux/r_cModelFluxErr > 20 "
        #f"AND i_cModelFlux/i_cModelFluxErr > 20 "
        #f"AND r_extendedness = 1"
    )
    query_lines.append("ORDER BY objectId ASC")

    query = "\n".join(query_lines)
    print(query)
    print()

    t0 = time.time()
    results= rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)




def select_highredshift_galaxies_grz_backup():

    title = 'dp1_ecdfs'

    query_lines = [
        "SELECT objectId, coord_ra, coord_dec, tract, patch, i_extendedness, "
]
    for filt in field_filters:
        query_lines.append(f"       {filt}_cModelFlux, {filt}_cModelFluxErr, ")
        query_lines.append(f"       {filt}_cModelMag AS {filt}mag, "
                       f"{filt}_cModelMagErr AS {filt}magErr, ")
    query_lines.append("       refExtendedness, refBand")

    query_lines.append("FROM dp1.Object")

    query_lines.append(
        f"WHERE CONTAINS(POINT('ICRS', coord_ra, coord_dec), "
        f"CIRCLE('ICRS', {ra_cen}, {dec_cen}, {radius})) = 1 "
        #f"AND r_cModelFlux/r_cModelFluxErr > 20 "
        #f"AND i_cModelFlux/i_cModelFluxErr > 20 "
        #f"AND r_extendedness = 1"
    )
    query_lines.append("ORDER BY objectId ASC")

    query = "\n".join(query_lines)
    print(query)
    print()

    t0 = time.time()
    results= rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)

    objtab = results.to_table()

    print()
    objtab.info(['attributes', 'stats'])

    objtab.write('objtab_dp1_ecdfs_galaxies.fits', overwrite=True)

    gmag = np.zeros(len(results['gmag']), dtype='float')
    S2N = np.asarray(results['g_cModelFlux']/results['g_cModelFluxErr'],
                     dtype='float')

    for i in range(len(gmag)):
        if S2N[i] < 1:
            gmag[i] = results['gmagErr'][i]
        else:
            gmag[i] = results['gmag'][i]

    GRlim = 1.5

    whlbg = np.where(((gmag-results['rmag']) > 1.5)
                     & ((results['gmag'] - results['rmag'])
                        > (results['rmag'] - results['zmag'] + 1.1))
                     & ((results['rmag'] - results['zmag']) < 1.0))[0]

    RmZ = np.arange(0.4, 1.1, 0.1)

    plt.plot(RmZ, RmZ + 1.1, '--', color='k')
    plt.plot([1., 1.], [2.1, 3], '--', color='k')
    plt.plot([-1, 0.4], [1.5, 1.5], '--', color='k')

    RminusZ = results['rmag'] - results['zmag']
    GminusR = gmag - results['rmag']

    label = 'all galaxies: ' + str(len(RminusZ))
    plt.plot(RminusZ, GminusR, '.', label=label, alpha=.1)

    label='color-selected z~4 candidates: ' + str(len(RminusZ[whlbg]))
    plt.plot(RminusZ[whlbg], GminusR[whlbg], 's', alpha=.5,
             label=label)

    plt.title(title)
    plt.xlabel('r - z')
    plt.ylabel('g - r')
    plt.xlim([-0.5, 3])
    plt.ylim([-0.5, 3])
    plt.legend()
    plt.show()

    return

def explore_objtab():


    return



def select_highredshift_quasars(mode='sync',
                                top=None,
                                allfields=False,
                                thintable=False):
    """
    fieldname?

    SELECT COUNT(*) FROM dp1.Object
    2299757 rows

    """

    thintable = True

    logger.info('')
    logging.info('')
    logger.info(f'Fieldname: {fieldname} top {top}')


    all = True
    allfields = True
    print('filters:', field_filters)

    query_lines = ["SELECT "]
    if top is not None:
        query_lines.append(f"       TOP {top} ")
    query_lines.append(f"       objectId, parentObjectId, ")
    query_lines.append(f"       coord_ra, coord_raErr, ")
    query_lines.append(f"       coord_dec, coord_decErr, ")
    query_lines.append(f"       tract, patch, ")
    query_lines.append(f"       refBand, refExtendedness, ")
    query_lines.append(f"       refSizeExtendedness, ")
    query_lines.append(f"       shape_flag, ")

    # note you need to put a variable with no comma at the end of
    # the SELECT series
    # also maybe a blank space at the end of each line
    for filt in field_filters:
        query_lines.append(f"       {filt}_ra, {filt}_raErr, ")
        query_lines.append(f"       {filt}_dec, {filt}_decErr, ")
        query_lines.append(f"       {filt}_i_flag, ")
        query_lines.append(f"       {filt}_extendedness, ")
        query_lines.append(f"       {filt}_extendedness_flag, ")
        if not thintable:
            query_lines.append(f"       {filt}_sizeExtendedness, ")
            query_lines.append(f"       {filt}_sizeExtendedness_flag, ")

        query_lines.append(f"       {filt}_psfFlux, ")
        query_lines.append(f"       {filt}_psfFluxErr, ")
        query_lines.append(f"       {filt}_psfFlux_flag, ")
        query_lines.append(f"       {filt}_free_psfFlux, ")
        query_lines.append(f"       {filt}_free_psfFluxErr, ")
        query_lines.append(f"       {filt}_free_psfFlux_flag, ")
        query_lines.append(f"       {filt}_psfMag, ")
        query_lines.append(f"       {filt}_psfMagErr, ")

        if not thintable:
            query_lines.append(f"       {filt}_cModelFlux, ")
            query_lines.append(f"       {filt}_cModelFlux_inner, ")
            query_lines.append(f"       {filt}_cModelFluxErr, ")
            query_lines.append(f"       {filt}_cModelMag, ")
            query_lines.append(f"       {filt}_cModelMagErr, ")
            query_lines.append(f"       {filt}_cModel_flag, ")
            query_lines.append(f"       {filt}_free_cModelFlux, ")
            query_lines.append(f"       {filt}_free_cModelFlux_inner, ")
            query_lines.append(f"       {filt}_free_cModelFluxErr, ")
            query_lines.append(f"       {filt}_free_cModelFlux_flag, ")

        query_lines.append(f"       {filt}_blendedness, ")
        query_lines.append(f"       {filt}_blendedness_flag, ")

    query_lines.append(f"       deblend_failed, deblend_masked, ")
    query_lines.append(f"       deblend_skipped ")
    query_lines.append("FROM dp1.Object ")

    if fieldname != 'All':
        query_lines.append(
            f"WHERE CONTAINS(POINT('ICRS', coord_ra, coord_dec), "
            f"CIRCLE('ICRS', {ra_centre}, {dec_centre}, {radius})) = 1 "
        )
    query_lines.append("ORDER BY objectId ASC")

    query = "\n".join(query_lines)

    print(query)
    print()

    t0 = time.time()
    print(f'Run query in mode: {mode} {fieldname} all: {all}')
    logger.info('')
    logging.info('')
    if mode == 'sync':
        result = rsp_tap.run_sync(query)

    if mode == 'async':
        result = run_async_query(service=rsp_tap, query=query)

    print(f'Number of rows returned {len(result)}')
    print('Elapsed time(secs): ',time.time() - t0)

    # objtab = result.to_table()

    print()
    result.info(['attributes'])
    # result.info(['attributes', 'stats'])

    # add some columns
    # refBand_S_N
    unique_refBands, counts = np.unique(result['refBand'],
                                        return_counts=True)

    print('Elapsed time(secs): ',time.time() - t0)
    print(f'Total number of rows: {len(result)}')
    result['refBand'].info(['attributes', 'stats'])
    for band, count in zip(unique_refBands, counts):
        print(f'Waveband count: {band}, {count}')
    logger.info('')

    # explore psfFlux
    refBand = result['refBand']
    print(type(result['u' + '_psfFlux']))
    nrows = len(result)
    refBand_psfFlux_S_N = np.zeros(nrows) * np.nan

    for irow, band in enumerate(refBand):
        refBand_psfFlux_S_N[irow] = (result[band + '_psfFlux'][irow] /
                             result[band + '_psfFluxErr'][irow])

    ndata_finite = np.count_nonzero(np.isfinite(refBand_psfFlux_S_N))
    print(f'S_N range: {np.min(refBand_psfFlux_S_N)} ' +
          f'{np.max(refBand_psfFlux_S_N)} ' +
          f'{np.median(refBand_psfFlux_S_N)} ' +
          f'{len(refBand_psfFlux_S_N)}')
    print(f'S_N range: {np.nanmin(refBand_psfFlux_S_N)} ' +
          f'{np.nanmax(refBand_psfFlux_S_N)} ' +
          f'{np.nanmedian(refBand_psfFlux_S_N)} ' +
          f'{len(refBand_psfFlux_S_N)}')
    print(f'Number of finite values: {ndata_finite}')

    print('Elapsed time(secs): ',time.time() - t0)

    u_psfFlux_S_N = result['u_psfFlux'] / result['u_psfFluxErr']
    g_psfFlux_S_N = result['g_psfFlux'] / result['g_psfFluxErr']
    r_psfFlux_S_N = result['r_psfFlux'] / result['r_psfFluxErr']
    i_psfFlux_S_N = result['i_psfFlux'] / result['i_psfFluxErr']
    z_psfFlux_S_N = result['z_psfFlux'] / result['z_psfFluxErr']
    y_psfFlux_S_N = result['y_psfFlux'] / result['y_psfFluxErr']

    result.add_column(refBand_psfFlux_S_N,
                       name='refBand_psfFlux_S_N',
                       index=10)

    result.add_columns([u_psfFlux_S_N, g_psfFlux_S_N,
                        r_psfFlux_S_N, i_psfFlux_S_N,
                        z_psfFlux_S_N, y_psfFlux_S_N],
                       names=['u_psfFlux_S_N', 'g_psfFlux_S_N',
                              'r_psfFlux_S_N', 'i_psfFlux_S_N',
                              'z_psfFlux_S_N', 'y_psfFlux_S_N'],
                       indexes=[12, 12, 12, 12, 12, 12])

    lu.plot_hist_psfFlux_S_N(table=result,
                             fluxtype='psfFlux',
                             fieldname=fieldname)


    fluxtype = 'free_psfFlux'
    refBand_free_psfFlux_S_N = np.zeros(nrows) * np.nan
    for irow, band in enumerate(refBand):
        refBand_free_psfFlux_S_N[irow] = (result[band + '_' + fluxtype][irow] /
                             result[band + '_' + fluxtype + 'Err'][irow])


    colname = 'refBand_free_psfFlux_S_N'
    ndata_finite = np.count_nonzero(np.isfinite(refBand_free_psfFlux_S_N))
    print(f'Number of finite values: {ndata_finite}')
    print(f'S_N range: {np.min(refBand_free_psfFlux_S_N)} ' +
          f'{np.max(refBand_free_psfFlux_S_N)} ' +
          f'{np.median(refBand_free_psfFlux_S_N)} ' +
          f'{len(refBand_free_psfFlux_S_N)}')
    print(f'S_N range: {np.nanmin(refBand_free_psfFlux_S_N)} ' +
          f'{np.nanmax(refBand_free_psfFlux_S_N)} ' +
          f'{np.nanmedian(refBand_free_psfFlux_S_N)} ' +
          f'{len(refBand_free_psfFlux_S_N)}')


    u_free_psfFlux_S_N = (result['u_' + fluxtype] /
                          result['u_' + fluxtype + 'Err'])
    g_free_psfFlux_S_N = (result['g_' + fluxtype] /
                          result['g_' + fluxtype + 'Err'])
    r_free_psfFlux_S_N = (result['r_' + fluxtype] /
                          result['r_' + fluxtype + 'Err'])
    i_free_psfFlux_S_N = (result['i_' + fluxtype] /
                          result['i_' + fluxtype + 'Err'])
    z_free_psfFlux_S_N = (result['z_' + fluxtype] /
                          result['z_' + fluxtype + 'Err'])
    y_free_psfFlux_S_N = (result['y_' + fluxtype] /
                          result['y_' + fluxtype + 'Err'])

    result.add_column(refBand_free_psfFlux_S_N,
                       name='refBand_' + fluxtype + '_S_N',
                       index=17)

    result.add_columns([u_free_psfFlux_S_N,
                        g_free_psfFlux_S_N,
                        r_free_psfFlux_S_N,
                        i_free_psfFlux_S_N,
                        z_free_psfFlux_S_N,
                        y_free_psfFlux_S_N],
                       names=['u_' + fluxtype + '_S_N',
                              'g_' + fluxtype + '_S_N',
                              'r_' + fluxtype + '_S_N',
                              'i_' + fluxtype + '_S_N',
                              'z_' + fluxtype + '_S_N',                                                       'y_' + fluxtype + '_S_N'],
                       indexes=[19, 19, 19, 19, 19, 19])

    logging.info('logging')
    logger.info('logger')
    help(result.info)
    print()
    lu.plot_hist_psfFlux_S_N(table=result,
                             fluxtype=fluxtype,
                             fieldname=fieldname)



    #

    result.info(['attributes'])
    outfile = 'DP1_' + fieldname + '_Object_tmp.fits'


    print(f'Saving: {outfile}')
    print('Elapsed time(secs): ',time.time() - t0)
    result.write(outfile, overwrite=True)
    print('Elapsed time(secs): ',time.time() - t0)

    return result

    run_async_query(service=rsp_tap)

    outfile = 'objtab_dp1_ecdfs_quasars.fits'
    print('Writing:', outfile)
    results.write('objtab_dp1_ecdfs_quasars.fits', overwrite=True)

    gmag = np.zeros(len(results['g_psfMag']), dtype='float')
    S2N = np.asarray(results['g_cModelFlux']/results['g_cModelFluxErr'],
                     dtype='float')

    for i in range(len(gmag)):
        if S2N[i] < 1:
            gmag[i] = results['g_psfMagErr'][i]
        else:
            gmag[i] = results['g_psfMag'][i]

    GRlim = 1.5

    whlbg = np.where(((gmag-results['r_psfMag']) > 1.5)
                     & ((results['g_psfMag'] - results['r_psfMag'])
                        > (results['r_psfMag'] - results['z_psfMag'] + 1.1))
                     & ((results['r_psfMag'] - results['z_psfMag']) < 1.0))[0]

    RmZ = np.arange(0.4, 1.1, 0.1)

    plt.plot(RmZ, RmZ + 1.1, '--', color='k')
    plt.plot([1., 1.], [2.1, 3], '--', color='k')
    plt.plot([-1, 0.4], [1.5, 1.5], '--', color='k')

    RminusZ = results['r_psfMag'] - results['z_psfMag']
    GminusR = gmag - results['r_psfMag']

    label = 'all galaxies: ' + str(len(RminusZ))
    plt.plot(RminusZ, GminusR, '.', label=label, alpha=.1)

    label='color-selected z~4 candidates: ' + str(len(RminusZ[whlbg]))
    plt.plot(RminusZ[whlbg], GminusR[whlbg], 's', alpha=.5,
             label=label)

    plt.title(title)
    plt.xlabel('r - z')
    plt.ylabel('g - r')
    plt.xlim([-0.5, 3])
    plt.ylim([-0.5, 3])
    plt.legend()
    plt.show()


    return


def mk_logger(prefix=None):

    logging.getLogger(__name__)

    logger_timestamp = time.strftime('%Y%m%dT%H%M', time.gmtime())

    if prefix is None:
        prefix = ''
    # set up file and screen logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s.%(msecs)03d %(name)-12s %(module)s %(funcName)s %(lineno)d %(levelname)-8s %(message)s',
        datefmt='%y-%m-%dT%H:%M:%S',
        filename=prefix + logger_timestamp + '.log',
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

    logging.info('Username: ' + username)
    logging.info(__name__)
    logging.info(__file__)

    logger  = logging

    return logger


if __name__ == "__main__":

    config = configparser.ConfigParser()

    timestamp = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
    filename_timestamp = time.strftime('%Y%m%dT%H%M', time.gmtime())

    username = getpass.getuser()
    print('__name__:', __name__)


    #from mk_logger import *
    logger = mk_logger()

    t0 = time.time()

    configfile = 'lsst_dp1.cfg'
    logger.info('Read configfile: ' + configfile)
    config.read(configfile)
    fieldnames = ['47Tuc', 'ECDFS', 'EDFS', 'LELF', 'FDSG']

    fieldname = fieldnames[1]
    sectionName = fieldname

    ntop = None
    # ntop = 100000
    # fieldname = 'All'

    logger.info(f'Fieldname: {fieldname}')
    radius = 1.0
    ra_centre = 53.2
    dec_centre = -28.1

    if fieldname !="All":

        field_ra_centre = config.get(sectionName, 'ra_centre')
        logger.info(f'Field RA centre {field_ra_centre}')

        field_dec_centre = config.get(sectionName, 'dec_centre')
        logger.info(f'Field Dec centre {field_dec_centre}')

        # target_ra = 53.2
        # target_dec = -28.1

        ra_centre = field_ra_centre
        dec_centre = field_dec_centre

        radius = 1.0


    plotfile_prefix = 'explore_dp1_' + fieldname

    run_xmatch_milliquas = True
    run_selfmatch = True
    UseCachedData = False

    t0 = time.time()
    rsp_tap = connect_rsp_tap()
    print('Elapsed time(secs): ',time.time() - t0, '\n')


    get_summary_info = False
    if get_summary_info:

        t0 = time.time()
        query_Object()
        print('Elapsed time(secs): ',time.time() - t0, '\n')

        t0 = time.time()
        query_CoaddPatches()
        print('Elapsed time(secs): ',time.time() - t0, '\n')

        t0 = time.time()
        query_schema(rsp_tap=rsp_tap)
        print('Elapsed time(secs): ',time.time() - t0, '\n')

        t0 = time.time()
        query_visits()
        print('Elapsed time(secs): ',time.time() - t0, '\n')


    # results= rsp_tap.run_sync(query)
    # query_ObsCore()
    # sys.exit()



    # select_highredshift_quasars(mode='async', top=10000)
    # query_object(mode='async', top=None)
    select_highredshift_quasars(mode='async', top=ntop)

    # sys.exit()

    run_ccdvisits = False
    if run_ccdvisits:
        t0 = time.time()
        query_ccdvisits(plotfile_prefix=plotfile_prefix)
        print('Elapsed time(secs): ',time.time() - t0, '\n')


    run_explore_objtab = False
    if run_explore_objtab:
        t0 = time.time()
        query_objtab()
        print('Elapsed time(secs): ',time.time() - t0, '\n')

    run_select_hzg = False
    if run_select_hzg:
        t0 = time.time()
        select_highredshift_galaxies_grz()
        print('Elapsed time(secs): ',time.time() - t0, '\n')

    sys.exit()

    """
    4. Object (detections)
    """


    # def explore_Magnitude_distributions()




    query = """SELECT objectId, coord_ra, coord_dec, ebv,
            {}_psfMag, {}_cModelMag, {}_psfMag, {}_cModelMag,
            {}_psfMag, {}_cModelMag, {}_psfMag, {}_cModelMag,
            {}_psfMag, {}_cModelMag, {}_psfMag, {}_cModelMag,
            refExtendedness
            FROM dp1.Object
            WHERE CONTAINS(POINT('ICRS', coord_ra, coord_dec),
                  CIRCLE('ICRS', {}, {}, {})) = 1
            ORDER BY objectId ASC
            """.format(field_filters[0], field_filters[0], field_filters[1], field_filters[1],
                       field_filters[2], field_filters[2], field_filters[3], field_filters[3],
                       field_filters[4], field_filters[4], field_filters[5], field_filters[5],
                       ra_cen, dec_cen, radius)
    print(query)

    t0 = time.time()
    results= rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)

    objtab = results.to_table()

    print()
    objtab.info(['attributes', 'stats'])

    objtab.write('objtab_dp1_ecdfs.fits')



    query = """
            SELECT visit, band, expMidptMJD, expMidpt
            FROM dp1.Visit
            WHERE CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {}, {}, {}))=1
            ORDER BY expMidptMJD
            """.format(ra_cen, dec_cen, radius)
    t0 = time.time()
    results = rsp_tap.run_sync(query)
    print('Elapsed time(secs): ',time.time() - t0)

    visits = results.to_table()

    visits.info(['attributes', 'stats'])

    print(f"Total number of visits: {len(visits)}")

    all_bands = np.array([visit['band'] for visit in visits], dtype=str)
    unique_bands, counts = np.unique(all_bands, return_counts=True)

    for band, count in zip(unique_bands, counts):
        print(band, count)

    # CcdVisits
    query = """
            SELECT visitId, ra, dec, band, magLim, seeing
            FROM dp1.CcdVisit
            WHERE CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {}, {}, {}))=1
            ORDER BY visitId
            """.format(ra_cen, dec_cen, radius)

    results = rsp_tap.run_sync(query)

    ccd_visits = results.to_table()


    for band in field_filters:
        print(f"Median in {band}-band: {np.ma.median(ccd_visits[ccd_visits['band'] == band]['seeing']): .2f}")


    min_seeing = np.min(ccd_visits['seeing'])
    max_seeing = np.max(ccd_visits['seeing'])
    seeing_bins = np.arange(min_seeing, max_seeing, 0.05)

    fig = plt.figure(figsize=(7, 5))
    for band in field_filters:
        print(f"Median in {band}-band: \
                {np.ma.median(ccd_visits[ccd_visits['band'] == band]['seeing']): .2f}")
        n, bins, patches = plt.hist(ccd_visits[ccd_visits['band'] == band]['seeing'],
                                    seeing_bins, histtype='step',
                                    linewidth=2,
                                    color=filter_colors[band],
                                    label=band)
        for patch in patches:
            patch.set_linestyle(filter_linestyles[band])

    plt.legend(loc='upper right')
    plt.xlabel('PSF FWHM (arcsec)')
    plt.ylabel('Number of visits')
    plt.minorticks_on()
    plt.show()


    min_mag = np.min(ccd_visits['magLim'])
    max_mag = np.max(ccd_visits['magLim'])
    maglim_bins = np.arange(min_mag, max_mag, 0.1)

    fig = plt.figure(figsize=(7, 5))

    for band in field_filters:
        print(f"Median in {band}-band: \
                {np.ma.median(ccd_visits[ccd_visits['band'] == band]['magLim']): .2f}")
        n, bins, patches = plt.hist(ccd_visits[ccd_visits['band'] == band]['magLim'],
                                    maglim_bins, histtype='step',
                                    linewidth=2,
                                    color=filter_colors[band],
                                    label=band)
        for patch in patches:
            patch.set_linestyle(filter_linestyles[band])

    plt.legend(loc='upper left')
    plt.xlabel('limiting magnitude')
    plt.ylabel('Number of visits')
    plt.minorticks_on()
    plt.show()

    """
    4. Object (detections)
    """




    # def explore_Magnitude_distributions()
