"""Utilities for the gene expression visualization portfolio project."""

from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version("gene-expression-visuals")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = ["__version__"]

