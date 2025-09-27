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

"""


import os
import sys
import time

import argparse
import logging

import numpy as np

import pyvo

from astropy.table import Table
import astropy.units as u


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Explore LSST RSP TAP service')

    parser.add_argument('-s', '--service',
                       choices=['US', 'UK', 'VSA', 'WSA'],
                       default='US',
                       help='Choose TAP service location (default: US)')

    return parser.parse_args()

def get_tap_service(service_name=None):
    """

    RSP US or UK
    VSA or WSA
    CDS ?
    your choice here

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

    print(f'Using {service_name} TAP Service URL: {tap_service_url}')

    tap_table_list = ['tap_schema']

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


def run_Query(service=None, query=None):

    print()
    print(f'query: {query}')
    print()

    result = service.run_sync(query)

    table = result.to_table()

    table.info(['attributes', 'stats'])

    return result

def make_query_region_count(table_name=None, adql_contain=False,
                            ra_centre=None, dec_centre=None,
                            radius=None, width=None, height=None):

    if adql_contain:
        query = """
            SELECT COUNT(*)
            FROM dp1.CoaddPatches
            WHERE CONTAINS(POINT('ICRS', s_ra, s_dec), CIRCLE('ICRS', {}, {}, {}))=1
            ORDER BY lsst_tract
            """.format(ra_centre, dec_centre, radius)

    query = f'SELECT COUNT(*) \n' + \
        f'FROM {table_name} \n' + \
        f'WHERE \n' + \
        f'    (ra >= 0.0) AND (ra <= +0.5) AND \n' + \
        f'    (dec >= 0.0) AND (dec <= +0.5)'

    query = f'SELECT COUNT(*) \n' + \
        f'FROM {table_name}'


    return query


def explore_vsa_wsa(service=None):
    """

    """
    print(service)
    help(service)

    t0 = time.time()

    table_names_vsa = ['VHSDR6.vhsSource', 'VIDEODR5.videoSource',
                       'VIKINGDR4.vikingSource']
                       # 'VMCDR5.vmcSource', 'VVVDR5.vvvSource']

    table_names_wsa = ['UKIDSSDR11PLUS.lasSource',
                       'UKIDSSDR11PLUS.dxsSource',
                       'UKIDSSDR11PLUS.gcsSource',
                       'UKIDSSDR11PLUS.gpsSource',
                       'UKIDSSDR11PLUS.udsSource',
                       'UHSDR2.uhsSource']

    for itable, table_name in enumerate(table_names_vsa):
        print(f'itable: {table_name}')
        query = make_query_region_count(table_name=table_name)

        result = run_Query(service=service, query=query)
        print(result)

        print(f'Elapsed time(secs): {time.time() - t0}\n')

    return


def main():

    t0 = time.time()

    # Parse command line arguments
    args = parse_arguments()
    print(f'args.service: {args.service}')

    service_name = args.service
    service = get_tap_service(service_name=service_name)
    print(f'Elapsed time(secs): {time.time() - t0}\n')


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
    print(f"Service URL: {service.baseurl}")
    print(f"Service type: {type(service)}")

    # Check if service has other URL-related attributes
    if hasattr(service, 'url'):
        print(f"URL attribute: {service.url}")

    input('Enter any key to continue... ')

    DEBUG = False
    if DEBUG:
        help(service)

    #query = "SELECT * FROM tap_schema."
    #result = run_Query(query=query)
    #print(result)

    query = "SELECT * FROM tap_schema.schemas"
    schemas = run_Query(service=service, query=query)
    print()
    input('Enter any key to continue... ')
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
        schema_name = 'VIKING'

    if args.service == 'VSA':
        schema_name = 'VIKINGDR4' # ['VHSDR6', 'VIDEODR5', 'VVVDR5']
        schema_name = 'VHSDR6' # ['VHSDR6', 'VIDEODR5', 'VVVDR5']

    if args.service == 'WSA':
        schema_name = 'UKIDSSDR11PLUS' # ['UHSDR2', 'UHSDR3']
        # schema_name = 'UHSDR2' # ['UHSDR2', 'UHSDR3']
        # schema_name = 'UHSDR3' # ['UHSDR2', 'UHSDR3']
        schema_names = ['UKIDSSDR11PLUS', 'UHSDR2']

    query = f"SELECT * FROM tap_schema.tables WHERE schema_name = '{schema_name}'"

    result = run_Query(service=service, query=query)
    print()
    for irow, row in enumerate(result):
        print(irow, row['table_name'])

    # sys.exit()

    print(type(service))
    print(service)
    help(service)
    table_name = 'dp1.Object'
    table_name = 'VHSDR6.vhsSource'


    query = f"SELECT column_name FROM tap_schema.columns WHERE table_name = '{table_name}'"

    result = run_Query(service=service, query=query)
    print(result)

    colnames_radec = ['ra', 'dec']
    query = make_query_region_count(table_name=table_name)
    print(f'query: \n{query}')
    result = run_Query(service=service, query=query)
    print(result)


    explore_vsa_wsa(service=service)

    sys.exit()

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
    main()
