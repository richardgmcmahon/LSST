#!/usr/bin/env python
"""
Example script showing how to get the number of rows and columns 
in a VizieR table using astroquery v0.4.12.dev255

This demonstrates using the VizierClass to query VizieR catalogs
and determine the dimensions of the returned table.

Usage:
    python vizier_table_dimensions_example.py
"""

from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord


def get_vizier_table_dimensions_basic():
    """
    Basic example: Query VizieR and get table dimensions.
    
    Returns the number of rows and columns in the result table.
    """
    print("\n=== Basic Example: Get VizieR Table Dimensions ===\n")
    
    # Create a VizieR query object
    v = Vizier()
    
    # Set row limit (optional, default is 50)
    v.ROW_LIMIT = 100
    
    # Query a catalog - example using 2MASS catalog
    # You can use catalog name or identifier
    catalog = "II/246"  # 2MASS Point Source Catalog
    
    # Define a position to search around
    coord = SkyCoord(ra=10.68458, dec=41.26917, unit=(u.deg, u.deg), frame='icrs')
    
    # Query VizieR around this position with a radius
    result = v.query_region(coord, radius=5*u.arcmin, catalog=catalog)
    
    # Check if results were returned
    if result:
        # result is a TableList - get the first table
        table = result[0]
        
        # Get number of rows (length of table)
        n_rows = len(table)
        
        # Get number of columns
        n_cols = len(table.colnames)
        
        print(f"Catalog: {catalog}")
        print(f"Number of rows: {n_rows}")
        print(f"Number of columns: {n_cols}")
        print(f"Column names: {table.colnames}")
        
        return n_rows, n_cols
    else:
        print("No results returned")
        return 0, 0


def get_vizier_table_dimensions_advanced():
    """
    Advanced example: Query VizieR with specific columns and get dimensions.
    
    Shows how to control which columns are returned and count them.
    """
    print("\n=== Advanced Example: Custom Columns ===\n")
    
    # Create VizieR query with specific columns
    v = Vizier(columns=['_RAJ2000', '_DEJ2000', 'Jmag', 'Hmag', 'Kmag'])
    v.ROW_LIMIT = 50
    
    # Query using catalog and position
    catalog = "II/246"
    coord = SkyCoord(ra=10.68458, dec=41.26917, unit=(u.deg, u.deg), frame='icrs')
    
    result = v.query_region(coord, radius=5*u.arcmin, catalog=catalog)
    
    if result:
        table = result[0]
        
        print(f"Catalog: {catalog}")
        print(f"Number of rows: {len(table)}")
        print(f"Number of columns: {len(table.colnames)}")
        print(f"Columns requested: {v.columns}")
        print(f"Actual columns returned: {table.colnames}")
        
        # You can also check table shape (rows, columns) if it's a numpy array
        print(f"\nTable info:")
        print(f"  Table type: {type(table)}")
        print(f"  Shape: ({len(table)} rows, {len(table.colnames)} columns)")
        
        return len(table), len(table.colnames)
    else:
        print("No results returned")
        return 0, 0


def get_vizier_table_dimensions_by_name():
    """
    Example: Query VizieR by catalog name and get dimensions.
    
    Different ways to query VizieR catalogs.
    """
    print("\n=== Query by Catalog Name ===\n")
    
    v = Vizier()
    v.ROW_LIMIT = 100
    
    # Query by object name instead of coordinates
    result = v.query_object("M31", catalog=["II/246"])
    
    if result:
        table = result[0]
        
        print(f"Object: M31")
        print(f"Number of rows: {len(table)}")
        print(f"Number of columns: {len(table.colnames)}")
        
        # Access individual table properties
        print(f"\nTable metadata:")
        print(f"  Table name: {table.meta.get('name', 'N/A')}")
        print(f"  Description: {table.meta.get('description', 'N/A')}")
        
        return len(table), len(table.colnames)
    else:
        print("No results returned")
        return 0, 0


def get_all_catalogs_info():
    """
    Example: Get information about multiple catalogs returned.
    
    VizieR can return multiple tables for a single query.
    """
    print("\n=== Multiple Catalogs Example ===\n")
    
    v = Vizier()
    v.ROW_LIMIT = 20
    
    # Query without specifying catalog - may return multiple tables
    coord = SkyCoord(ra=10.68458, dec=41.26917, unit=(u.deg, u.deg), frame='icrs')
    result = v.query_region(coord, radius=5*u.arcmin, catalog="II/246")
    
    print(f"Number of tables returned: {len(result)}")
    
    for i, table in enumerate(result):
        print(f"\nTable {i+1}:")
        print(f"  Number of rows: {len(table)}")
        print(f"  Number of columns: {len(table.colnames)}")
        print(f"  Sample columns: {table.colnames[:5]}")


def main():
    """
    Main function demonstrating different ways to get VizieR table dimensions.
    """
    print("=" * 70)
    print("VizieR Table Dimensions with astroquery.vizier.VizierClass")
    print("=" * 70)
    
    # Run all examples
    try:
        get_vizier_table_dimensions_basic()
        get_vizier_table_dimensions_advanced()
        get_vizier_table_dimensions_by_name()
        get_all_catalogs_info()
        
        print("\n" + "=" * 70)
        print("Summary: How to get rows and columns")
        print("=" * 70)
        print("""
Key Methods:
1. len(table)           - Returns number of rows
2. len(table.colnames)  - Returns number of columns
3. table.colnames       - Returns list of column names
4. table.info()         - Shows detailed table information

Example Code:
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    
    v = Vizier()
    v.ROW_LIMIT = 100
    coord = SkyCoord(ra=10.68, dec=41.27, unit='deg')
    result = v.query_region(coord, radius=5*u.arcmin, catalog="II/246")
    
    if result:
        table = result[0]
        n_rows = len(table)
        n_cols = len(table.colnames)
        print(f"Rows: {n_rows}, Columns: {n_cols}")
        """)
        
    except Exception as e:
        print(f"\nError running examples: {e}")
        print("Make sure astroquery is installed: pip install astroquery")


if __name__ == "__main__":
    main()
