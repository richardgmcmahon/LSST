"""
Module for creating test astropy tables for testing purposes.
"""

import numpy as np
from astropy.table import Table


def create_test_table(n_rows=100, 
                     column_names=None,
                     column_ranges=None,
                     seed=None):
    """
    Create a test astropy table with random data.
    
    Parameters
    ----------
    n_rows : int, optional
        Number of rows in the table (default: 100)
    column_names : list of str, optional
        Names of columns to create. Default is ['redshift', 'mag_g', 'mag_r']
    column_ranges : list of tuple, optional
        Ranges for each column as (min, max). If a single tuple is provided,
        it will be used for all columns. Default is (-10, 10) for all columns.
    seed : int, optional
        Random seed for reproducibility
    
    Returns
    -------
    astropy.table.Table
        Test table with random floating point data
    
    Examples
    --------
    >>> # Create a basic test table with defaults
    >>> table = create_test_table()
    
    >>> # Create a table with specific columns
    >>> table = create_test_table(
    ...     n_rows=100,
    ...     column_names=['redshift', 'mag_g', 'mag_r'],
    ...     column_ranges=[(-10, 10)]
    ... )
    
    >>> # Create a table with different ranges for each column
    >>> table = create_test_table(
    ...     column_names=['z', 'flux_g', 'flux_r'],
    ...     column_ranges=[(0, 3), (0, 1000), (0, 1000)]
    ... )
    """
    
    # Set default column names if not provided
    if column_names is None:
        column_names = ['redshift', 'mag_g', 'mag_r']
    
    # Set default column ranges if not provided
    if column_ranges is None:
        column_ranges = [(-10, 10)]
    
    # If only one range is provided, use it for all columns
    if len(column_ranges) == 1:
        column_ranges = column_ranges * len(column_names)
    
    # Validate that we have the same number of ranges as columns
    if len(column_ranges) != len(column_names):
        raise ValueError(
            f"Number of column_ranges ({len(column_ranges)}) must match "
            f"number of column_names ({len(column_names)}) or be 1"
        )
    
    # Set random seed if provided
    if seed is not None:
        np.random.seed(seed)
    
    # Create the table
    table = Table()
    
    # Add columns with random data in the specified ranges
    for col_name, (min_val, max_val) in zip(column_names, column_ranges):
        # Generate random floating point data in the range [min_val, max_val]
        data = np.random.uniform(min_val, max_val, size=n_rows)
        table[col_name] = data
    
    return table


if __name__ == '__main__':
    # Example usage
    print("Creating test table with default parameters...")
    table = create_test_table()
    print(table)
    print(f"\nTable info:")
    print(f"  Number of rows: {len(table)}")
    print(f"  Column names: {table.colnames}")
    print(f"  Column dtypes: {[table[col].dtype for col in table.colnames]}")
    
    print("\n" + "="*60)
    print("Creating test table with custom parameters...")
    table2 = create_test_table(
        n_rows=50,
        column_names=['redshift', 'mag_g', 'mag_r'],
        column_ranges=[(0, 5), (15, 25), (15, 25)],
        seed=42
    )
    print(table2[:10])  # Show first 10 rows
    print(f"\nTable info:")
    print(f"  Number of rows: {len(table2)}")
    print(f"  Column ranges:")
    for col in table2.colnames:
        print(f"    {col}: [{table2[col].min():.2f}, {table2[col].max():.2f}]")
