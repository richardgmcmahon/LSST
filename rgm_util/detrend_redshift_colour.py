"""
Module for detrending redshift-colour relationships.
"""

import logging

from .get_footnote import get_footnote

logger = logging.getLogger(__name__)


def detrend_redshift_colour(table=None, 
                            redshift_col='z',
                            colour_col='colour',
                            detrend_method='ndimage',
                            verbose=False,
                            **kwargs):
    """
    Detrend redshift-colour relationships.
    
    This function removes systematic trends in colour as a function of redshift,
    which can be useful for identifying outliers or unusual objects.
    
    Parameters
    ----------
    table : astropy.table.Table
        Input table containing the data
    redshift_col : str, optional
        Name of the redshift column (default: 'z')
    colour_col : str, optional
        Name of the colour column (default: 'colour')
    detrend_method : str, optional
        Method to use for detrending. Options: 'ndimage', 'polynomial' (default: 'ndimage')
    verbose : bool, optional
        Enable verbose output (default: False)
    **kwargs : dict
        Additional keyword arguments passed to the detrending function
    
    Returns
    -------
    table : astropy.table.Table
        Table with added detrended colour column
    """
    
    # Import numpy and scipy here to avoid import errors if not installed
    try:
        import numpy as np
        from scipy import ndimage
    except ImportError as e:
        logger.error(f"Required dependency not available: {e}")
        raise
    
    if table is None:
        raise ValueError("Input table cannot be None")
    
    if redshift_col not in table.colnames:
        raise ValueError(f"Column '{redshift_col}' not found in table")
    
    if colour_col not in table.colnames:
        raise ValueError(f"Column '{colour_col}' not found in table")
    
    logger.info(f"Detrend using scipy {detrend_method}")
    
    # Get the data
    redshift = np.array(table[redshift_col])
    colour = np.array(table[colour_col])
    
    # Remove NaN values for processing
    mask = np.isfinite(redshift) & np.isfinite(colour)
    
    if detrend_method == 'ndimage':
        # Use scipy ndimage for smoothing
        # This is a placeholder implementation
        # Real implementation would depend on specific requirements
        from scipy.ndimage import median_filter
        
        # Sort by redshift
        sort_idx = np.argsort(redshift[mask])
        sorted_z = redshift[mask][sort_idx]
        sorted_colour = colour[mask][sort_idx]
        
        # Apply median filter
        window_size = kwargs.get('window_size', 51)
        trend = median_filter(sorted_colour, size=window_size, mode='nearest')
        
        # Compute detrended values
        detrended = sorted_colour - trend
        
        # Map back to original order
        detrended_colour = np.full_like(colour, np.nan)
        detrended_colour[mask] = detrended[np.argsort(sort_idx)]
        
    elif detrend_method == 'polynomial':
        # Polynomial fitting
        degree = kwargs.get('degree', 2)
        coeffs = np.polyfit(redshift[mask], colour[mask], degree)
        trend = np.polyval(coeffs, redshift)
        detrended_colour = colour - trend
        
    else:
        raise ValueError(f"Unknown detrend method: '{detrend_method}'")
    
    # Add detrended column to table
    detrended_col_name = f'{colour_col}_detrended'
    table[detrended_col_name] = detrended_colour
    
    if verbose:
        logger.info(f"Added column '{detrended_col_name}' to table")
    
    footnote_text = get_footnote()
    if verbose:
        logger.debug(f"Generated footnote: {footnote_text}")
    
    return table
