"""Build and package the Windows executable in the runtime checkout."""

import argparse
import hashlib
import importlib.metadata
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


def collect_notices(project):
    import imageio_ffmpeg
    destination = project / 'build' / 'third_party'
    destination.mkdir(parents=True, exist_ok=True)
    (destination / 'PROJECT-LICENSE.txt').write_text(
        'Voronoi Cell ODE — https://github.com/LBK888/Voronoi_Cell_ODE\n'
        'The project README identifies the application source as CC BY-NC 4.0.\n'
        'https://creativecommons.org/licenses/by-nc/4.0/\n'
        'Third-party components retain their own licenses, included separately.\n', encoding='utf-8')
    packages = []
    for distribution in importlib.metadata.distributions():
        name = distribution.metadata['Name']
        packages.append(f'{name}=={distribution.version}')
        for entry in distribution.files or []:
            parts = Path(str(entry)).parts
            if not any(part.endswith('.dist-info') for part in parts):
                continue
            if not ('licenses' in parts or Path(str(entry)).name.upper().startswith(('LICENSE', 'COPYING', 'NOTICE'))):
                continue
            source = Path(distribution.locate_file(entry))
            if source.is_file():
                target = destination / name / Path(*parts[1:])
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
    (destination / 'BUILD-ENVIRONMENT.txt').write_text('\n'.join(sorted(packages)) + '\n', encoding='utf-8')
    encoder = imageio_ffmpeg.get_ffmpeg_exe()
    for option, name in (('-version', 'FFMPEG-VERSION.txt'), ('-L', 'FFMPEG-LICENSE.txt')):
        completed = subprocess.run([encoder, option], capture_output=True, text=True, check=True)
        (destination / name).write_text(completed.stdout + completed.stderr, encoding='utf-8')
    (destination / 'FFMPEG-SOURCES.txt').write_text(
        'The bundled FFmpeg is a separate executable supplied by imageio-ffmpeg.\n'
        'Its build details and license are in FFMPEG-VERSION.txt and FFMPEG-LICENSE.txt.\n'
        'Binary distributor, build scripts and corresponding source links:\n'
        'https://www.gyan.dev/ffmpeg/builds/\n'
        'https://github.com/GyanD/codexffmpeg\n'
        'https://ffmpeg.org/download.html\n'
        'Python wrapper / binary selection: https://github.com/imageio/imageio-ffmpeg\n',
        encoding='utf-8',
    )
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', default='2026.11.08')
    args = parser.parse_args()
    if sys.platform != 'win32':
        parser.error('Run this build with 64-bit Windows Python')
    project = Path(__file__).resolve().parents[1]
    notices = collect_notices(project)
    # Prevent unrelated Conda/document-tool DLLs on the user's PATH from
    # replacing Windows system libraries during dependency collection.
    environment = os.environ.copy()
    windows = Path(os.environ.get('SystemRoot', 'C:/Windows'))
    environment['PATH'] = os.pathsep.join(map(str, [
        Path(sys.executable).parent, Path(sys.base_prefix),
        Path(sys.base_prefix) / 'DLLs',
        Path(sys.prefix) / 'Lib/site-packages/PyQt6/Qt6/bin',
        windows / 'System32', windows,
    ]))
    environment.pop('PYTHONPATH', None)
    subprocess.run([sys.executable, '-m', 'PyInstaller', '--noconfirm', '--clean',
                    '--distpath', str(project / 'dist'),
                    '--workpath', str(project / 'build' / 'pyinstaller'),
                    str(project / 'packaging' / 'windows.spec')], cwd=project, env=environment, check=True)
    executable = project / 'dist' / 'VoronoiCellODE.exe'
    smoke = project / 'build' / 'frozen-smoke'
    # The executable must work without Python, FFmpeg or build tools on PATH.
    environment['PATH'] = os.pathsep.join(map(str, [windows / 'System32', windows]))
    subprocess.run([str(executable), '--smoke-test', str(smoke)], cwd=project,
                   env=environment, check=True, timeout=180)
    report = smoke / 'smoke-report.json'
    if not report.is_file():
        raise RuntimeError('Frozen application did not produce a smoke-test report')
    import json
    if json.loads(report.read_text(encoding='utf-8')).get('status') != 'passed':
        raise RuntimeError('Frozen application smoke test failed')
    portable = project / 'dist' / f'VoronoiCellODE-{args.version}-windows-x64.zip'
    with zipfile.ZipFile(portable, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        archive.write(executable, 'VoronoiCellODE.exe')
        archive.write(project / 'docs' / 'windows-release.md', 'README.md')
        for folder in (notices, project / 'docs', project / 'skills'):
            prefix = 'third_party' if folder == notices else folder.name
            for item in folder.rglob('*'):
                if item.is_file() and '__pycache__' not in item.parts:
                    archive.write(item, str(Path(prefix) / item.relative_to(folder)))
    checksums = project / 'dist' / 'SHA256SUMS.txt'
    lines = []
    for path in (executable, portable):
        with path.open('rb') as handle:
            lines.append(f'{hashlib.file_digest(handle, "sha256").hexdigest()}  {path.name}\n')
    checksums.write_text(''.join(lines), encoding='utf-8')
    print(f'Built and smoke-tested: {executable}\nPortable ZIP: {portable}\nChecksums: {checksums}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
