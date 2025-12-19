"""
Module for generating footnotes for plots and figures.
"""

import os
import datetime


def get_footnote(timestamp=None, filename=None, additional_info=None):
    """
    Generate a footnote string for plots and figures.
    
    Parameters
    ----------
    timestamp : str, optional
        Timestamp to include in the footnote. If None, current time is used.
    filename : str, optional
        Filename to include in the footnote. If None, calling file name is used.
    additional_info : str, optional
        Additional information to include in the footnote.
    
    Returns
    -------
    str
        Formatted footnote string
    """
    if timestamp is None:
        timestamp = datetime.datetime.now().strftime('%y-%m-%dT%H:%M:%S.%f')[:-3]
    
    parts = []
    if timestamp:
        parts.append(timestamp)
    if filename:
        parts.append(os.path.basename(filename))
    if additional_info:
        parts.append(additional_info)
    
    return ': '.join(parts) if parts else ''
