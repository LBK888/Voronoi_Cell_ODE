"""Separate bundled resources from persistent, writable user settings."""

import os
from pathlib import Path
import sys


def resource_path(name):
    root = Path(getattr(sys, '_MEIPASS', Path(__file__).resolve().parents[1]))
    return root / name


def data_path(name):
    override = os.environ.get('VORONOI_CELL_ODE_DATA_DIR')
    if override:
        root = Path(override)
    elif getattr(sys, 'frozen', False):
        root = Path(os.environ.get('LOCALAPPDATA') or Path.home() / 'AppData' / 'Local') / 'VoronoiCellODE'
    else:
        return resource_path(name)
    root.mkdir(parents=True, exist_ok=True)
    return root / name


def configure_ffmpeg():
    """Use the bundled encoder when running without a system FFmpeg install."""
    import matplotlib
    bundled = resource_path('bin/ffmpeg.exe')
    if bundled.is_file():
        matplotlib.rcParams['animation.ffmpeg_path'] = str(bundled)
