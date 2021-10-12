try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version("healsparse")
except PackageNotFoundError:
    # package is not installed
    pass

from .fdest import FgcmDesTransmission
