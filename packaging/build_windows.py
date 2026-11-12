"""Build and package the Windows executable in the runtime checkout."""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import mkdtemp
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
    parser.add_argument('--version', default='2026.11.12')
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
    smoke = Path(mkdtemp(prefix='frozen-smoke-', dir=project / 'build'))
    install = smoke / 'portable'
    install.mkdir()
    installed_exe = install / executable.name
    try:
        os.link(executable, installed_exe)
    except OSError:
        shutil.copy2(executable, installed_exe)
    shutil.copytree(project / 'templates', install / 'templates')
    # The executable must work without Python, FFmpeg or build tools on PATH.
    environment['PATH'] = os.pathsep.join(map(str, [windows / 'System32', windows]))
    validate_package(project, smoke, installed_exe, environment)
    portable = project / 'dist' / f'VoronoiCellODE-{args.version}-windows-x64.zip'
    with zipfile.ZipFile(portable, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        archive.write(executable, 'VoronoiCellODE.exe')
        archive.write(project / 'docs' / 'windows-release.md', 'README.md')
        for folder in (notices, project / 'docs', project / 'skills', project / 'templates'):
            prefix = 'third_party' if folder == notices else folder.name
            for item in folder.rglob('*'):
                if item.is_file() and '__pycache__' not in item.parts:
                    archive.write(item, str(Path(prefix) / item.relative_to(folder)))
    templates_zip = project / 'dist' / 'VoronoiCellODE-templates.zip'
    with zipfile.ZipFile(templates_zip, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted((project / 'templates').glob('*.json')):
            archive.write(path, f'templates/{path.name}')
    checksums = project / 'dist' / 'SHA256SUMS.txt'
    lines = []
    for path in (executable, portable, templates_zip):
        with path.open('rb') as handle:
            lines.append(f'{hashlib.file_digest(handle, "sha256").hexdigest()}  {path.name}\n')
    checksums.write_text(''.join(lines), encoding='utf-8')
    print(f'Built and smoke-tested: {executable}\nPortable ZIP: {portable}\n'
          f'Template ZIP: {templates_zip}\nSmoke reports: {smoke}\nChecksums: {checksums}')
    return 0


def validate_package(project, smoke, executable, environment):
    model = project / 'skills/paper-to-voronoi-ode/assets/hill-demo.json'
    subprocess.run([str(executable), '--smoke-test', str(smoke), '--smoke-model', str(model)],
                   cwd=smoke, env=environment, check=True, timeout=180)
    report = json.loads((smoke / 'smoke-report.json').read_text(encoding='utf-8'))
    expected_templates = {path.name for path in (project / 'templates').glob('*.json')}
    expected_templates.add('Smoke Template.json')
    if report.get('status') != 'passed' or {Path(path).name for path in report['templates']} != expected_templates:
        raise RuntimeError('Frozen GUI did not discover the external and saved templates.')

    environment = {**environment, 'VORONOI_CELL_ODE_DATA_DIR': str(smoke / 'settings')}
    sweep = smoke / 'sweep.json'
    sweep.write_text(json.dumps({
        'x': {'target': 'alpha', 'minimum': 1, 'maximum': 2, 'count': 2},
        'y': {'kind': 'initial', 'target': 1, 'minimum': 0.1, 'maximum': 0.5, 'count': 2},
    }), encoding='utf-8')
    for command, args in (
        ('list-templates', []),
        ('simulate', [report['saved_template'], '--output', str(smoke / 'simulation')]),
        ('experiment', [report['saved_template'], '--output', str(smoke / 'sweep'), '--settings', str(sweep)]),
    ):
        response = smoke / f'{command}-report.json'
        subprocess.run([str(executable), '--cli', '--report', str(response), command, *args],
                       cwd=smoke, env=environment, check=True, timeout=180)
        result = json.loads(response.read_text(encoding='utf-8'))
        if result.get('status') != 'ok' or result.get('failed_runs', 0):
            raise RuntimeError(f'Frozen CLI {command} failed: {result}')
        if command == 'list-templates':
            if {Path(item['path']).name for item in result['templates']} != expected_templates:
                raise RuntimeError('Frozen CLI did not discover the external templates.')
        elif not all(Path(path).is_file() for path in result['files']):
            raise RuntimeError(f'Frozen CLI {command} did not produce its output files.')


if __name__ == '__main__':
    raise SystemExit(main())
