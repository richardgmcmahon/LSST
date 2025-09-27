"""

https://www.ivoa.net/documents/TAP/

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


def run_Query(service=None, query=None):

    print()
    print(f'query: {query}')
    print()

    result = service.run_sync(query)

    table = result.to_table()

    table.info(['attributes', 'stats'])

    return result


def main():

    # Parse command line arguments
    args = parse_arguments()
    print(f'args.service: {args.service}')

    # Set service URL based on command line argument
    if args.service == 'VSA':
        TAP_SERVICE = 'http://tap.roe.ac.uk/vsa/'
    elif args.service == 'WSA':
        TAP_SERVICE = 'http://tap.roe.ac.uk/wsa/'
    elif args.service == 'US':
        TAP_SERVICE = 'https://data.lsst.cloud/api/tap'
    elif args.service == 'UK':
        TAP_SERVICE = 'https://rsp.lsst.ac.uk/api/tap'
    else:
        TAP_SERVICE = None

    print(f'Using {args.service} TAP Service URL: {TAP_SERVICE}')

    tap_table_list = ['tap_schema']

    homedir = os.path.expanduser('~')
    token_file = None
    if args.service == 'US':
        token_file = os.path.join(homedir,'.rsp-tap.token')
    if args.service == 'UK':
        token_file = os.path.join(homedir,'.rsp_uk-tap.token')

    rsp_tap = pyvo.dal.TAPService(TAP_SERVICE)
    if token_file is not None:
        with open(token_file, 'r') as f:
            token_str = f.readline()

        cred = pyvo.auth.CredentialStore()
        cred.set_password("x-oauth-basic", token_str)
        credential = cred.get("ivo://ivoa.net/sso#BasicAA")
        rsp_tap = pyvo.dal.TAPService(TAP_SERVICE, session=credential)

    service = rsp_tap

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

    if args.service == 'WSA':
        schema_name = 'UKIDSSDR11PLUS' # ['UHSDR2', 'UHSDR3']
        schema_name = 'UHSDR2' # ['UHSDR2', 'UHSDR3']
        schema_name = 'UHSDR3' # ['UHSDR2', 'UHSDR3']

    query = f"SELECT * FROM tap_schema.tables WHERE schema_name = '{schema_name}'"

    result = run_Query(service=service, query=query)
    print()
    for irow, row in enumerate(result):
        print(irow, row['table_name'])

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
