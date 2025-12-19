"""
Module for writing RA/Dec coordinates to CSV files.
"""

import logging
import os

logger = logging.getLogger(__name__)


def write_radec_csvfile(table=None, filename='radec.csv', 
                       ra_col='ra', dec_col='dec', **kwargs):
    """
    Write RA and Dec coordinates to a CSV file.
    
    Parameters
    ----------
    table : astropy.table.Table
        Input table with coordinates
    filename : str, optional
        Output filename
    ra_col : str, optional
        Name of RA column
    dec_col : str, optional
        Name of Dec column
    **kwargs : dict
        Additional keyword arguments
    
    Returns
    -------
    filename : str
        Path to the created file
    """
    logger.info(f"Writing RA/Dec to {filename}")
    
    if table is None:
        raise ValueError("Input table cannot be None")
    
    # Write a simple CSV file
    with open(filename, 'w') as f:
        f.write(f"{ra_col},{dec_col}\n")
        for row in table:
            f.write(f"{row[ra_col]},{row[dec_col]}\n")
    
    logger.info(f"Wrote {len(table)} coordinates to {filename}")
    
    return filename
