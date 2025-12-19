"""
Module for getting git hash information.
"""

import logging
import subprocess

logger = logging.getLogger(__name__)


def get_githash(short_hash=True):
    """
    Get the current git commit hash.
    
    Parameters
    ----------
    short_hash : bool, optional
        If True, return short hash (7 chars), otherwise full hash
    
    Returns
    -------
    str
        Git commit hash or 'unknown' if not in a git repository
    """
    try:
        if short_hash:
            result = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                                  capture_output=True, text=True, check=True, timeout=5)
        else:
            result = subprocess.run(['git', 'rev-parse', 'HEAD'],
                                  capture_output=True, text=True, check=True, timeout=5)
        
        git_hash = result.stdout.strip()
        logger.debug(f"Git hash: {git_hash}")
        return git_hash
        
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        logger.warning("Could not determine git hash")
        return 'unknown'
