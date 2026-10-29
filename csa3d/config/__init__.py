"""Configuration package for CSA-3D.

Handles data paths, environment variables, and data loading.
"""

from .settings import get_config, set_data_dir, get_data_dir, validate_config

__all__ = [
    "get_config",
    "set_data_dir",
    "get_data_dir",
    "validate_config",
]
