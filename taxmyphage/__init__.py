"""
This module is used to set the version of the package.
It uses the version information from the package metadata.
"""

from importlib.metadata import version

__version__ = version(__package__)
