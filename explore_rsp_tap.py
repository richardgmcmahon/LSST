"""

https://www.ivoa.net/documents/TAP/

https://wiki.ivoa.net/twiki/bin/view/IVOA/ADQL

https://data.lsst.cloud/

https://sdm-schemas.lsst.io/

https://sdm-schemas.lsst.io/dp1.html#Object

https://dp0-2.lsst.io/data-access-analysis-tools/api-intro.html#use-of-pyvo-with-the-rsp-tap-service

https://community.lsst.org/t/will-there-be-external-tap-access-to-rsp-dp0-2-tables/6660

https://rsp.lsst.ac.uk/nb/user/richardgmcmahon/lab/tree/notebooks/rsp-uk-notebooks/VISTA-HSC_tutorials/5_SimpleTAPDemo.ipynb


NIR Fusion relevant
WSA: UKIDSSDR11PLUS, UHSDR2/UHSDR3
VSA:


UHSDR3: 559,972,658, 1,258,772,895
UHSDR2: 513,163,485, 820,052,738
UHSDR1: 480,796,968, 481,056,015


IRSA
https://irsa.ipac.caltech.edu/docs/program_interface/astropy_TAP.html


ESO:
https://archive.eso.org/programmatic/

https://archive.eso.org/programmatic/HOWTO/jupyter/how_to_query_for_reduced_data/ESO_How_to_query_for_reduced_data.html


CDS:
http://tapvizier.u-strasbg.fr/tapvizier/tap/
https://cdsarc.cds.unistra.fr/vizier/notebook.gml?source=B/eso

"""

import os
import sys
import time

import argparse
import logging

import inspect
import time
from datetime import datetime

from SimpleLogger import SimpleLogger

import numpy as np

import pyvo

from astropy.table import Table
import astropy.units as u



class SimpleLogger:
    """Simple logger with automatic line number and function name detection."""

    """
    logger written by Claude with this prompt

    can you add a logger that writes to the screen the line number and
    the function name and add it before the function call main()

    e.g. logger.info('Hello World')

    Simple Logger with File Support

    Console Logging:
        logger.info('Hello World')
        logger.debug('Debug message')
        logger.warning('Warning message')
        logger.error('Error message')

    File Logging:
        logger.file('Info to file and console')
        logger.file('Info to file only', also_console=False)
        logger.file('Custom file', filename='custom.log')
        logger.file_debug('Debug to file only')  # No console by default
        logger.file_error('Error to file and console')

    File Management:
        logger.set_log_file('new_default.log')
        logger.clear_log_file()  # Clear current log file
        logger.clear_log_file('specific.log')  # Clear specific file


    # ========================================================================
    # LOGGER USAGE EXAMPLES
    # ========================================================================

    # Basic console logging
    logger.info('Hello World')
    logger.debug('This is a debug message')
    logger.warning('This is a warning')
    logger.error('This is an error message')

    # File logging examples
    logger.file('This goes to file and console')
    logger.file('This goes to file only', also_console=False)
    logger.file('Custom log file', filename='custom.log')

    # Debug and error file logging
    logger.file_debug('Debug message to file only')  # No console by default
    logger.file_error('Error message to file and console')

    # File management
    logger.set_log_file('my_custom.log')
    logger.file('Now using custom log file')

    # Different logger configurations
    print("\n" + "="*60)
    print("DIFFERENT LOGGER CONFIGURATIONS:")
    print("="*60)

    # Quiet logger (no timestamps or levels)
    quiet_logger = SimpleLogger(show_time=False, show_level=False, log_file='quiet.log')
    quiet_logger.info('Clean output without timestamps or levels')
    quiet_logger.file('Clean file logging too')

    # Debug logger (with all details)
    debug_logger = SimpleLogger(show_time=True, show_level=True, log_file='debug.log')
    debug_logger.debug('Detailed debug information')
    debug_logger.file_debug('Debug info to file')

    print("\n" + "="*60)
    print("MAIN FUNCTION EXECUTION:")
    print("="*60)

    # Main function call with logging
    logger.info('About to call main() function')
    logger.file('Starting main() execution')

    try:
        main()
    except Exception as e:
        logger.error(f'main() failed: {e}')
        logger.file_error(f'main() execution failed: {e}')

    logger.info('Finished calling main() function')
    logger.file('main() execution completed')

    print("\n" + "="*60)
    print("ADDITIONAL USAGE EXAMPLES:")
    print("="*60)

    # Example with variables
    service_name = 'VSA'
    table_count = 42
    query_time = 1.234

    logger.info(f'Processing {service_name} service')
    logger.info(f'Found {table_count} tables in {query_time:.3f} seconds')
    logger.file(f'Performance: {table_count} tables processed in {query_time} seconds')

    # Example with error handling
    try:
        # Simulate some operation
        result = 10 / 0  # This will cause an error
    except ZeroDivisionError as e:
        logger.error(f'Division error occurred: {e}')
        logger.file_error(f'Mathematical error: {e}')

    # Example with step-by-step logging
    logger.info('Starting multi-step process')
    logger.debug('Step 1: Initialize variables')
    logger.debug('Step 2: Connect to service')
    logger.debug('Step 3: Execute queries')
    logger.info('Multi-step process completed')

    print("\n" + "="*60)
    print("Check the following log files:")
    print("- tap_exploration.log (main log)")
    print("- custom.log (custom file example)")
    print("- quiet.log (quiet logger example)")
    print("- debug.log (debug logger example)")
    print("="*60)

    """

    import os
    import sys
    import time

    import logging
    import inspect
    from datetime import datetime

    def __init__(self, show_time=True, show_level=True, log_file='script.log'):
        """Initialize the logger.

        Parameters
        ----------
        show_time : bool, optional
            Whether to include timestamp in output
        show_level : bool, optional
            Whether to include log level in output
        log_file : str, optional
            Default filename for file logging
        """
        self.show_time = show_time
        self.show_level = show_level
        self.log_file = log_file

    def _log(self, level, message, to_file=False, filename=None, also_console=True):
        """Internal logging method."""
        # Get the calling frame (skip this _log method)
        frame = inspect.currentframe().f_back.f_back

        # Extract information
        line_number = frame.f_lineno
        function_name = frame.f_code.co_name
        filename_code = frame.f_code.co_filename.split('/')[-1]  # Just filename, not full path

        # Format timestamp
        timestamp = ""
        if self.show_time:
            timestamp = f"[{datetime.now().strftime('%H:%M:%S')}] "

        # Format level
        level_str = ""
        if self.show_level:
            level_str = f"{level.upper()}: "

        # Format the log message
        log_msg = f"{timestamp}{level_str}{filename_code}:{line_number} in {function_name}() - {message}"

        if to_file:
            # Write to file
            file_to_use = filename if filename else self.log_file
            try:
                with open(file_to_use, 'a', encoding='utf-8') as f:
                    f.write(log_msg + '\n')
            except Exception as e:
                print(f"ERROR: Could not write to log file {file_to_use}: {e}")

            # Also print to console if requested
            if also_console:
                print(log_msg)
        else:
            # Just print to console
            print(log_msg)

    def info(self, message):
        """Log an info message."""
        self._log("info", message)

    def debug(self, message):
        """Log a debug message."""
        self._log("debug", message)

    def error(self, message):
        """Log an error message."""
        self._log("error", message)

    def warning(self, message):
        """Log a warning message."""
        self._log("warning", message)

    def warn(self, message):
        """Log a warning message (alias for warning)."""
        self.warning(message)

    def file(self, message, filename=None, also_console=True):
        """Log an info message to file.

        Parameters
        ----------
        message : str
            Message to log
        filename : str, optional
            Specific filename to write to. If None, uses default log_file
        also_console : bool, optional
            Whether to also print to console. Default is True
        """
        self._log("info", message, to_file=True, filename=filename, also_console=also_console)

    def file_debug(self, message, filename=None, also_console=False):
        """Log a debug message to file (console=False by default for debug)."""
        self._log("debug", message, to_file=True, filename=filename, also_console=also_console)

    def file_error(self, message, filename=None, also_console=True):
        """Log an error message to file."""
        self._log("error", message, to_file=True, filename=filename, also_console=also_console)

    def set_log_file(self, filename):
        """Set the default log file name."""
        self.log_file = filename

    def clear_log_file(self, filename=None):
        """Clear/truncate the log file."""
        file_to_clear = filename if filename else self.log_file
        try:
            with open(file_to_clear, 'w', encoding='utf-8') as f:
                f.write(f"Log cleared at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            print(f"Log file {file_to_clear} cleared")
        except Exception as e:
            print(f"ERROR: Could not clear log file {file_to_clear}: {e}")



def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Explore LSST RSP TAP service')

    parser.add_argument('-s', '--service',
                       choices=['US', 'UK', 'VSA', 'WSA', 'IRSA', 'ESO','CDS'],
                       default='US',
                       help='Choose TAP service location (default: US)')

    return parser.parse_args()


def test_cds():
    """
    https://www.euro-vo.org/the-cds-tutorial/
    """

    tap_vizier = pyvo.dal.TAPService('http://tapvizier.u-strasbg.fr/TAPVizieR/tap')
    query = """SELECT  *  FROM tap_schema.tables
               WHERE table_name LIKE '%arp%' """

    catalog_list = tap_vizier.search(query).to_table()

    print(catalog_list['table_name', 'description'])

    query = """SELECT  *  FROM tap_schema.tables
               WHERE table_name LIKE '%eso_arc%' """

    catalog_list = tap_vizier.search(query).to_table()

    help(catalog_list)
    catalog_list.info()
    print(catalog_list)
    print(catalog_list['table_name', 'description'])



    sys.exit()

    tap_url = "http://simbad.cds.unistra.fr/simbad/sim-tap"
    service = pyvo.dal.TAPService(tap_url)
    simbad = service
    print(service)
    help(service)
    print([tab_name for tab_name in service.tables.keys()])

    print()

    tap_url = 'https://tapvizier.u-strasbg.fr/adql'
    tap_url = 'http://tapvizier.u-strasbg.fr/TAPVizieR/tap'
    service = pyvo.dal.TAPService(tap_url)
    vizier = service
    print(service)

    #for key, item in service.tables.items():
    #    print(key, item)

    print()
    for key in service.tables.keys():
        print(key)

    # print([tab_name for tab_name in service.tables.keys()])

    print()

    return



def get_tap_service(service_name=None):
    """

    RSP US or UK
    VSA or WSA
    CDS ?
    your choice here

    """

    # Dictionary mapping service names to service URLs
    service_urls = {
        'US': 'https://data.lsst.cloud/api/tap',
        'UK': 'https://rsp.lsst.ac.uk/api/tap',
        'VSA': 'http://tap.roe.ac.uk/vsa/',
        'WSA': 'http://tap.roe.ac.uk/wsa/',
        'ESO': 'http://archive.eso.org/tap_obs',
        'IRSA': 'https://irsa.ipac.caltech.edu/TAP',
        'CDS': 'http://tapvizier.u-strasbg.fr/TAPVizieR/tap/',
        'NED': 'https://ned.ipac.caltech.edu/tap/',
    }

    # Set service URL based on command line argument
    if service_name in service_urls:
        tap_service_url = service_urls[service_name]
    else:
        raise ValueError(f"Unknown service: {service_name}. Available services: {list(service_urls.keys())}")

    print(f'tap service_url: {tap_service_url}')

    # sys.exit()

    """
    # Set service URL based on command line argument
    tap_service_url = None
    if service_name == 'VSA':
        tap_service_url = 'http://tap.roe.ac.uk/vsa/'
    elif service_name == 'WSA':
        tap_service_url = 'http://tap.roe.ac.uk/wsa/'
    elif service_name == 'US':
        tap_service_url = 'https://data.lsst.cloud/api/tap'
    elif service_name == 'UK':
        tap_service_url = 'https://rsp.lsst.ac.uk/api/tap'
    else:
        tap_service_url = None
    """

    print(f'Using {service_name} TAP Service URL: {tap_service_url}')


    homedir = os.path.expanduser('~')
    token_file = None
    if service_name == 'US':
        token_file = os.path.join(homedir,'.rsp-tap.token')
    if service_name == 'UK':
        token_file = os.path.join(homedir,'.rsp_uk-tap.token')

    rsp_tap = pyvo.dal.TAPService(tap_service_url)
    if token_file is not None:
        with open(token_file, 'r') as f:
            token_str = f.readline()

        cred = pyvo.auth.CredentialStore()
        cred.set_password("x-oauth-basic", token_str)
        credential = cred.get("ivo://ivoa.net/sso#BasicAA")
        rsp_tap = pyvo.dal.TAPService(tap_service_url, session=credential)

    service = rsp_tap

    return service


def get_service_info(service=None):
    """


    """
    t0 = time.time()

    # Handle services that do not expose limit information
    try:
        maxrec = service.maxrec
        print(f'service.maxrec: {maxrec}')
    except Exception as e:
        print(f'service.maxrec: Not available ({e})')

    try:
        hardlimit = service.hardlimit
        print(f'service.hardlimit: {hardlimit}')
    except Exception as e:
        print(f'service.hardlimit: Not available ({e})')
    print(f'Elapsed time(secs): {time.time() - t0}\n')

    # Print service info
    # this is quite slow
    get_service_info = False
    if get_service_info:
        print(f'Get the service description using service.describe')
        print(f'{service.describe(width=80)}')
        print(f'Elapsed time(secs): {time.time() - t0}\n')


    print(f'Get the service examples')
    print(f'{service.examples}')
    print(f'Elapsed time(secs): {time.time() - t0}\n')


    get_service_info = True
    if get_service_info:
        try:
            print(f"Service info: {service.info}")
        except:
            print("No service info available")
    print(f'Elapsed time(secs): {time.time() - t0}\n')

    # Print service capabilities (if available)
    try:
        print(f"Service description: {service.description}")
    except:
        print("No service description available")
    print(f'Elapsed time(secs): {time.time() - t0}\n')

    service_baseurl = service.baseurl
    print(f"service.baseurl: {service.baseurl}")
    print(f"Service type: {type(service)}")

    # Check if service has other URL-related attributes
    if hasattr(service, 'url'):
        print(f"URL attribute: {service.url}")


    return



def run_Query(service=None, query=None):

    print()
    print(f'query: {query}')
    print()

    result = service.run_sync(query)

    table = result.to_table()

    table.info(['attributes', 'stats'])

    return result


def make_query_region(table_name=None,
                      colnames_radec=None,
                      count=False,
                      count_all=False,
                      adql_contain=False,
                      primary=True,
                      ra_centre=180.0, dec_centre=0.0,
                      radius=1.0,
                      width=None, height=None):
    """Construct an ADQL spatial query for astronomical catalogs.

    Creates ADQL queries for spatial searches in astronomical survey tables,
    supporting both geometric containment functions and simple coordinate
    range queries.

    For VISTA and WFCAM primary sources are defined by
     (priOrSec<=0 OR priOrSec=frameSetID)

    Parameters
    ----------
    table_name : str, optional
        Name of the table to query (e.g., 'dp1.CoaddPatches', 'VHSDR6.vhsSource').
        If None, a default table may be used.
    count : bool, optional
        If True, return COUNT(*) query instead of selecting all columns.
        Default is False.
    adql_contain : bool, optional
        If True, use ADQL CONTAINS function for spatial matching.
        If False, use simple coordinate range constraints.
        Default is False.
    ra_centre : float, optional
        Right ascension of search center in decimal degrees (ICRS).
        Required when adql_contain=True.
    dec_centre : float, optional
        Declination of search center in decimal degrees (ICRS).
        Required when adql_contain=True.
    radius : float, optional
        Search radius in decimal degrees.
        Required when adql_contain=True.
    width : float, optional
        Width of rectangular search region in decimal degrees.
        Used when adql_contain=False.
    height : float, optional
        Height of rectangular search region in decimal degrees.
        Used when adql_contain=False.

    Returns
    -------
    str
        ADQL query string ready for execution via TAP service.

    Examples
    --------
    Create a circular search query using CONTAINS:

    >>> query = make_query_region(
    ...     table_name='dp1.CoaddPatches',
    ...     adql_contain=True,
    ...     ra_centre=150.0,
    ...     dec_centre=2.0,
    ...     radius=0.1
    ... )

    Create a simple count query:

    >>> query = make_query_region(
    ...     table_name='VHSDR6.vhsSource',
    ...     count=True
    ... )

    Notes
    -----
    The function currently has incomplete logic and will always return
    a simple COUNT(*) query regardless of parameters. This appears to
    be a work-in-progress implementation.

    ADQL spatial functions follow IVOA standards:
    - POINT('ICRS', ra, dec) creates a point in ICRS coordinates
    - CIRCLE('ICRS', ra, dec, radius) creates a circular region
    - CONTAINS(geometry1, geometry2) tests spatial containment
    """

    import astropy.units as u
    from astropy.coordinates import SkyCoord

    if width is None:
        width = 2.0 * radius
    if height is None:
        height = width

    if adql_contain:
        query = """
            SELECT COUNT(*)
            FROM dp1.CoaddPatches
            WHERE CONTAINS(POINT('ICRS', s_ra, s_dec), CIRCLE('ICRS', {}, {}, {}))=1
            ORDER BY lsst_tract
            """.format(ra_centre, dec_centre, radius)

    DEBUG = True

    if not adql_contain and count_all:
                query = f'SELECT COUNT(*) \n' + \
            f'FROM {table_name}'

    if not adql_contain and not count_all:
        dec_limits = \
                [dec_centre - (0.5 * height),
                 dec_centre + (0.5 * height)]
        maxabsdec = np.max(np.abs(dec_limits))
        if DEBUG:
            print(f'maxabsdec: {maxabsdec}')

        cosdec_correction = 1.0/(np.cos(np.deg2rad(maxabsdec)))
        if DEBUG:
            print(f'cosdec_correction: {cosdec_correction}')

        ra_limits = \
            [ra_centre - (0.5 * width * cosdec_correction),
             ra_centre + (0.5 * width * cosdec_correction)]
        print(f'ra_limits: {ra_limits}')
        print(f'ra_limits: {dec_limits}')

        c1 = SkyCoord(ra_limits[0]*u.deg, dec_limits[0]*u.deg,
                      frame='icrs')
        c2 = SkyCoord(ra_limits[1]*u.deg, dec_limits[0]*u.deg,
                      frame='icrs')
        sep = c1.separation(c2)
        print(f'RA range: {sep.degree} degrees')

        c1 = SkyCoord(ra_limits[0]*u.deg, dec_limits[1]*u.deg,
                      frame='icrs')
        c2 = SkyCoord(ra_limits[1]*u.deg, dec_limits[1]*u.deg,
                      frame='icrs')
        sep = c1.separation(c2)
        print(f'RA range: {sep.degree} degrees')

        query = f'SELECT COUNT(*) \n' + \
            f'FROM {table_name} \n' + \
            f'WHERE \n' + \
            f'    (ra >= {ra_limits[0]}) AND (ra <= {ra_limits[1]}) AND \n' + \
            f'    (dec >= {dec_limits[0]}) AND (dec <= {dec_limits[1]})'

    print(f'TAP ADQL query: \n {query}')

    return query


def explore_rsp(service=None):

    print(type(service))
    print(service)
    # help(service)

    t0 = time.time()

    if 'cloud' in service.baseurl.lower():
        table_names = ['dp1.CcdVisit',
                       'dp1.CoaddPatches',
                       'dp1.DiaObject',
                       'dp1.DiaSource',
                       'dp1.ForcedSource',
                       'dp1.ForcedSourceOnDiaObject',
                       'dp1.MPCORB',
                       'dp1.Object',
                       'dp1.Source',
                       'dp1.SSObject',
                       'dp1.SSSource',
                       'dp1.Visit']



    for itable, table_name in enumerate(table_names):
        print(f'itable: {table_name}')
        query = f'SELECT COUNT(*) FROM {table_name} '

        result = run_Query(service=service, query=query)
        print(result)
        print(f'Elapsed time(secs): {time.time() - t0}\n')

    return


def explore_vsa_wsa(service=None):
    """

    """
    print(type(service))
    print(service)
    # help(service)

    t0 = time.time()

    if 'vsa' in service.baseurl.lower():
        table_names = ['VHSDR6.vhsSource',
                       'VIDEODR5.videoSource',
                       'VIKINGDR4.vikingSource']
                       #'VMCDR5.vmcSource']
                       #'VVVDR4.vvvSource']

    if 'wsa' in service.baseurl.lower():
        table_names = ['UKIDSSDR11PLUS.lasSource',
                       'UKIDSSDR11PLUS.dxsSource',
                       'UKIDSSDR11PLUS.gcsSource',
                       'UKIDSSDR11PLUS.gpsSource',
                       'UKIDSSDR11PLUS.udsSource',
                       'UHSDR2.uhsSource']

    for itable, table_name in enumerate(table_names):
        print(f'itable: {table_name}')
        if count:
            query = f'SELECT COUNT(*) FROM {table_name})'

        if not count:
            query = make_query_region(
                table_name=table_name,
                adql_contain=False,
                count=True)



        result = run_Query(service=service, query=query)
        print(result)
        print(f'Elapsed time(secs): {time.time() - t0}\n')



    return table_names


def main():
    """


    """


    t0 = time.time()

    # Parse command line arguments
    logger.info('Parse command line arguments')
    args = parse_arguments()
    print(f'args.service: {args.service}')

    logger.info(f'Get the TAP service: {args.service}')
    service_name = args.service
    service = get_tap_service(service_name=service_name)
    print(type(service))
    print(service)
    # help(service)
    print(f'Elapsed time(secs): {time.time() - t0}\n')

    logger.info('Get some information about the TAP service')
    get_service_info(service=service)

    DEBUG = False
    if DEBUG:
        help(service)

    logger.info('\n')
    input('Enter any key to continue... ')


    #query = "SELECT * FROM tap_schema."
    #result = run_Query(query=query)
    #print(result)

    # tap_table_list = ['tap_schema']

    #this does not work
    #query = "SELECT * FROM tap_schema"
    #tap_schema = run_Query(service=service, query=query)
    #table_count = len(tap_schema)
    #logger.info(f'Number of rows in tap_schema in: {table_count}\n')
    #input('Enter any key to continue... ')

    logger.info('List the TAP service table schema')
    query = "SELECT * FROM tap_schema.schemas"
    schemas = run_Query(service=service, query=query)
    table_count = len(schemas)
    logger.info(f'Number of table schemas in TAP service: {table_count}\n')
    input('Enter any key to continue... ')

    logger.info(f'Print names of table schema names in TAP service')
    print(f'schema_name')
    for irow, row in enumerate(schemas):
        print(irow, row['schema_name'])

    #input('Enter any key to continue... ')
    #for irow, row in enumerate(schemas):
    #    print()
    #    print(irow, 'schema_name:', row['schema_name'])
    #    print('utype:', row['utype'])
    #    print('description:', row['description'])
    #    print('schema_index:', row['schema_index'])

    print()
    input('Enter any key to continue... ')

    query = "SELECT * FROM tap_schema.tables"
    tables = run_Query(service=service, query=query)
    print()
    for itable, table in enumerate(tables):
        print(itable, table['table_type'], table['table_name'])
        print(itable, table['description'])
        print()

    schema_name = None
    if args.service == 'US':
        schema_name = 'DP1'

    if args.service == 'UK':
        # schema_name = 'VIKING'
        schema_name = 'VIDEO'

    if args.service == 'VSA':
        # schema_name = 'VIKINGDR4' # ['VHSDR6', 'VIDEODR5', 'VVVDR5']
        schema_name = 'VHSDR6' # ['VHSDR6', 'VIDEODR5', 'VVVDR5']

    if args.service == 'WSA':
        schema_name = 'UKIDSSDR11PLUS' # ['UHSDR2', 'UHSDR3']
        # schema_name = 'UHSDR2' # ['UHSDR2', 'UHSDR3']
        # schema_name = 'UHSDR3' # ['UHSDR2', 'UHSDR3']
        schema_names = ['UKIDSSDR11PLUS', 'UHSDR2']

    logger.info(f'Explore schema name: {schema_name}')
    query = f"SELECT * FROM tap_schema.tables WHERE schema_name = '{schema_name}'"

    result = run_Query(service=service, query=query)
    print()
    for irow, row in enumerate(result):
        print(irow, row['table_name'])

    # sys.exit()


    # Example data table_name
    table_name = None
    if 'lsst.cloud' in service.baseurl.lower():
        table_name = 'dp1.Object'
    if 'lsst.ac.uk' in service.baseurl.lower():
        table_name = 'video.merged'
    if 'vsa' in service.baseurl.lower():
        table_name = 'VHSDR6.vhsSource'
    if 'wsa' in service.baseurl.lower():
        table_name = 'UKIDSSDR11PLUS.lasSource'
    if 'eso' in service.baseurl.lower():
        table_name = 'phase3v2.files'
    if 'eso' in service.baseurl.lower():
        table_name = 'B/eso/eso_arc'

    query = f"SELECT column_name FROM tap_schema.columns WHERE table_name = '{table_name}'"

    logger.info(f'run Query against TAP service table: {table_name}')
    result = run_Query(service=service, query=query)
    print(result)

    logger.info(f'run RA, Dec region query against TAP service table: {table_name}')
    colnames_radec = ['ra', 'dec']
    query = make_query_region(table_name=table_name, count_all=True,
                              colnames_radec=None)
    print(f'query: \n{query}')

    result = run_Query(service=service, query=query)
    print(result)

    print(f'service.baseurl: {service.baseurl.lower()}')
    if any(test in service.baseurl.lower() for test in ['vsa', 'wsa']):
        logger.info('Run explore_vsa_wsa(service=service)')
        explore_vsa_wsa(service=service)

    if any(test in service.baseurl.lower() for test in ['cloud']):
        logger.info('Run explore_rsp(service=service)')
        explore_rsp(service=service)


    input('Enter any key to continue... ')


    # sys.exit()

    for itable, table in enumerate(tables):
        print(itable, table['table_type'], table['table_name'])
        print(itable, table['description'])
        print()


    print('List just DP1 tables')
    input('Enter any key to continue... ')
    for itable, table in enumerate(tables):
        if ('dp1' in table['table_name']):
            print(itable, table['table_type'], table['table_name'])
            print(itable, table['description'])
            print()
    print()

    print('List just IVOA tables')
    input('Enter any key to continue... ')
    for itable, table in enumerate(tables):
        if ('ivoa' in table['table_name']):
            print(itable, table['table_type'], table['table_name'])
            print(itable, table['description'])
            print()
    print()

    print('List just tap_schema tables')
    input('Enter any key to continue... ')
    for itable, table in enumerate(tables):
        if ('tap_schema' in table['table_name']):
            print(itable, table['table_type'], table['table_name'])
            print(itable, table['description'])
            print()
    print()


    input('Enter any key to continue... ')

    query = "SELECT DISTINCT table_name  FROM tap_schema.columns"
    tables = run_Query(service=service, query=query)

    query = "SELECT * FROM tap_schema.columns"
    columns = run_Query(service=service, query=query)
    print()

    input('Enter any key to continue... ')
    for icol, column in enumerate(columns):
        print()
        print(icol, column['table_name'], column['column_name'])
        print('Description:', column['description'].strip(), ':')
        print('utype:', column['utype'].strip())
        print('ucd:', column['ucd'].strip())
        print('unit:',column['unit'].strip())
        print('datatype:', column['datatype'].strip())


    query = "SELECT * FROM tap_schema.keys"
    run_Query(service=service, query=query)

    query = "SELECT * FROM tap_schema.key_columns"
    run_Query(service=service, query=query)

    return

if __name__ == "__main__":

    # Create a global logger instance
    logger = SimpleLogger()

    # Alternative: create logger without timestamps for cleaner output
    # logger = SimpleLogger(show_time=False)

    # Alternative: create logger without level names
    # logger = SimpleLogger(show_level=False)


    #test_cds()
    #sys.exit()

    logger.info('About to call main() function\n')
    main()
