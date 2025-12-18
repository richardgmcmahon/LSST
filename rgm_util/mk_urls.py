"""
Module for creating URLs for various astronomical services.
"""

import logging

logger = logging.getLogger(__name__)


def mk_legacy_survey_url(ra, dec, layer='ls-dr9', zoom=14, size=512):
    """
    Create a Legacy Survey viewer URL.
    
    Parameters
    ----------
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    layer : str, optional
        Image layer to display
    zoom : int, optional
        Zoom level
    size : int, optional
        Image size in pixels
    
    Returns
    -------
    str
        URL to Legacy Survey viewer
    """
    url = f"https://www.legacysurvey.org/viewer?ra={ra}&dec={dec}&layer={layer}&zoom={zoom}"
    return url


def mk_sdss_url(ra, dec, scale=0.4):
    """
    Create an SDSS Navigate Tool URL.
    
    Parameters
    ----------
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    scale : float, optional
        Image scale in arcsec/pixel
    
    Returns
    -------
    str
        URL to SDSS Navigate Tool
    """
    url = f"http://skyserver.sdss.org/dr16/en/tools/chart/navi.aspx?ra={ra}&dec={dec}&scale={scale}"
    return url


def mk_aladin_url(ra, dec, fov=0.2):
    """
    Create an Aladin Lite URL.
    
    Parameters
    ----------
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    fov : float, optional
        Field of view in degrees
    
    Returns
    -------
    str
        URL to Aladin Lite
    """
    url = f"https://aladin.u-strasbg.fr/AladinLite/?target={ra}+{dec}&fov={fov}"
    return url


__all__ = [
    'mk_legacy_survey_url',
    'mk_sdss_url',
    'mk_aladin_url',
]
