# Build with: python -m PyInstaller packaging/windows.spec --noconfirm
from pathlib import Path
import imageio_ffmpeg

project = Path(SPECPATH).parent
notice_dir = project / 'build' / 'third_party'

analysis = Analysis(
    [str(project / 'windows_launcher.py')],
    pathex=[str(project)],
    binaries=[(imageio_ffmpeg.get_ffmpeg_exe(), 'bin')],
    datas=[
        (str(project / 'config_example.json'), '.'),
        (str(project / 'templates'), 'templates'),
        (str(project / 'skills'), 'skills'),
        (str(notice_dir), 'third_party'),
    ],
    hiddenimports=['openpyxl', 'h5py', 'PyQt6.Qsci', 'gui.node_editor',
                   'gui.code_editor_window', 'matplotlib.backends.backend_pdf',
                   'matplotlib.backends.backend_svg'],
    hookspath=[],
    hooksconfig={'matplotlib': {'backends': ['QtAgg', 'Agg']}},
    runtime_hooks=[],
    excludes=['tkinter', 'PyQt5', 'PySide2', 'PySide6', 'IPython', 'pytest'],
    noarchive=False,
)
# Give the bundled encoder a stable name used by gui.app_paths.
analysis.binaries = [
    ('bin/ffmpeg.exe', source, kind) if Path(name).parent.name == 'bin' and Path(name).name.startswith('ffmpeg-')
    else (name, source, kind)
    for name, source, kind in analysis.binaries
]
pyz = PYZ(analysis.pure)
exe = EXE(
    pyz, analysis.scripts, analysis.binaries, analysis.datas, [],
    name='VoronoiCellODE', debug=False, bootloader_ignore_signals=False,
    strip=False, upx=False, console=False, disable_windowed_traceback=False,
    argv_emulation=False, target_arch=None, codesign_identity=None,
    entitlements_file=None,
)
