"""
This module is used to set the version of the package.
It uses the version information from the package metadata.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version(__package__)
except PackageNotFoundError:
    __version__ = 'dev version'
