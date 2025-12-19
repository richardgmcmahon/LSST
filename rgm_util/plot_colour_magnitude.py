"""
Module for creating colour-magnitude plots.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def plot_colour_magnitude(table=None, mag_col='mag', colour_col='colour',
                          showplots=True, **kwargs):
    """
    Create a colour-magnitude diagram.
    
    Parameters
    ----------
    table : astropy.table.Table
        Input table containing photometric data
    mag_col : str, optional
        Name of magnitude column
    colour_col : str, optional
        Name of colour column
    showplots : bool, optional
        Whether to display the plot
    **kwargs : dict
        Additional keyword arguments for plotting
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    logger.info("Creating colour-magnitude diagram")
    
    if table is None:
        raise ValueError("Input table cannot be None")
    
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10, 8)))
    
    if mag_col in table.colnames and colour_col in table.colnames:
        ax.scatter(table[colour_col], table[mag_col], alpha=0.5, s=1)
        ax.set_xlabel(colour_col)
        ax.set_ylabel(mag_col)
        ax.invert_yaxis()
        ax.grid(True, alpha=0.3)
    
    if showplots:
        plt.show()
    
    return fig
