"""
Module for plotting RA/Dec coordinates on the sky.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def plot_radec(table=None, ra_col='ra', dec_col='dec', 
              showplots=True, **kwargs):
    """
    Plot RA and Dec coordinates on the sky.
    
    Parameters
    ----------
    table : astropy.table.Table
        Input table with coordinates
    ra_col : str, optional
        Name of RA column
    dec_col : str, optional
        Name of Dec column
    showplots : bool, optional
        Whether to display the plot
    **kwargs : dict
        Additional keyword arguments for plotting
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    logger.info("Plotting RA/Dec coordinates")
    
    if table is None:
        raise ValueError("Input table cannot be None")
    
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (12, 8)),
                          subplot_kw={'projection': kwargs.get('projection', None)})
    
    if ra_col in table.colnames and dec_col in table.colnames:
        ax.scatter(table[ra_col], table[dec_col], alpha=0.5, s=1)
        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')
        ax.grid(True, alpha=0.3)
        ax.invert_xaxis()  # RA increases to the left
    
    if showplots:
        plt.show()
    
    return fig
