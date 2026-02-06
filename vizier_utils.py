"""
Utility functions for working with VizieR tables from astroquery.

This module provides helper functions for working with VizieR catalog queries
using astroquery v0.4.12+ and the VizierClass.
"""

from astropy.table import Table


def get_vizier_table_dimensions(table):
    """
    Get the number of rows and columns in a VizieR table.
    
    This is a helper function for working with astroquery VizieR results.
    When using astroquery.vizier.Vizier to query catalogs, the result is
    typically a TableList containing astropy Table objects.
    
    Parameters
    ----------
    table : astropy.table.Table
        An astropy Table object returned from a VizieR query
        
    Returns
    -------
    tuple of (int, int)
        Number of rows and number of columns (n_rows, n_cols)
        
    Examples
    --------
    >>> from astroquery.vizier import Vizier
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> from vizier_utils import get_vizier_table_dimensions
    >>> 
    >>> # Query VizieR
    >>> v = Vizier()
    >>> v.ROW_LIMIT = 100
    >>> coord = SkyCoord(ra=10.68, dec=41.27, unit='deg')
    >>> result = v.query_region(coord, radius=5*u.arcmin, catalog="II/246")
    >>> 
    >>> if result:
    >>>     table = result[0]
    >>>     n_rows, n_cols = get_vizier_table_dimensions(table)
    >>>     print(f"Table has {n_rows} rows and {n_cols} columns")
    
    Notes
    -----
    This function works with any astropy Table object, not just VizieR results.
    
    For VizieR tables, you can also directly access:
    - len(table) - number of rows
    - len(table.colnames) - number of columns
    - table.colnames - list of column names
    - table.info() - detailed table information
    
    See Also
    --------
    vizier_table_dimensions_example.py : Complete examples of VizieR queries
    """
    if not hasattr(table, '__len__') or not hasattr(table, 'colnames'):
        raise TypeError("Input must be an astropy Table or compatible object")
    
    n_rows = len(table)
    n_cols = len(table.colnames)
    
    return n_rows, n_cols


def print_vizier_table_info(table, name=None):
    """
    Print comprehensive information about a VizieR table.
    
    Parameters
    ----------
    table : astropy.table.Table
        An astropy Table object
    name : str, optional
        Name to display for the table
        
    Examples
    --------
    >>> from astroquery.vizier import Vizier
    >>> from vizier_utils import print_vizier_table_info
    >>> v = Vizier()
    >>> result = v.query_object("M31", catalog=["II/246"])
    >>> if result:
    >>>     print_vizier_table_info(result[0], name="2MASS M31")
    """
    n_rows, n_cols = get_vizier_table_dimensions(table)
    
    if name:
        print(f"\n=== Table: {name} ===")
    else:
        print(f"\n=== Table Information ===")
    
    print(f"Number of rows: {n_rows}")
    print(f"Number of columns: {n_cols}")
    print(f"Column names: {table.colnames}")
    
    # Show metadata if available
    if hasattr(table, 'meta') and table.meta:
        print(f"\nMetadata:")
        for key, value in table.meta.items():
            if key in ['name', 'description', 'ID']:
                print(f"  {key}: {value}")
