"""
rgm_util - Utility functions for astronomical data processing

Note: Some functions require additional dependencies (numpy, scipy, astropy, matplotlib).
Import individual modules as needed to avoid dependency issues.
"""

# Lazy imports - only import what's actually used
# This prevents import errors if dependencies like numpy/scipy aren't installed

__all__ = [
    'detrend_redshift_colour',
    'get_footnote',
    'explore_flux_fluxerr',
    'plot_colour_magnitude',
    'xmatch_ls',
    'write_radec_csvfile',
    'get_githash',
    'plot_radec',
    'create_test_table',
]


def __getattr__(name):
    """Lazy import of submodules."""
    if name == 'detrend_redshift_colour':
        from .detrend_redshift_colour import detrend_redshift_colour
        return detrend_redshift_colour
    elif name == 'get_footnote':
        from .get_footnote import get_footnote
        return get_footnote
    elif name == 'explore_flux_fluxerr':
        from .explore_flux_fluxerr import explore_flux_fluxerr
        return explore_flux_fluxerr
    elif name == 'plot_colour_magnitude':
        from .plot_colour_magnitude import plot_colour_magnitude
        return plot_colour_magnitude
    elif name == 'xmatch_ls':
        from .xmatch_ls import xmatch_ls
        return xmatch_ls
    elif name == 'write_radec_csvfile':
        from .write_radec_csvfile import write_radec_csvfile
        return write_radec_csvfile
    elif name == 'get_githash':
        from .get_githash import get_githash
        return get_githash
    elif name == 'plot_radec':
        from .plot_radec import plot_radec
        return plot_radec
    elif name == 'create_test_table':
        from .create_test_table import create_test_table
        return create_test_table
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
