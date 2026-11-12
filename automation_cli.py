"""Run local model jobs in an independent process with JSON reports."""

import argparse
from contextlib import redirect_stdout
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from gui.app_paths import data_path
from gui.model_files import model_grid, model_request, read_model, template_files, write_json


def cli_enabled():
    try:
        settings = json.loads(data_path('automation.json').read_text(encoding='utf-8'))
        return isinstance(settings, dict) and settings.get('cli_enabled') is True
    except (OSError, ValueError):
        return False


def set_cli_enabled(enabled):
    write_json(data_path('automation.json'), {'cli_enabled': enabled})


class CliParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


def cli_parser():
    parser = CliParser(prog='VoronoiCellODE --cli', description='Run local simulation jobs.')
    parser.add_argument('--report', type=Path, help='Write the JSON response to this file.')
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('list-templates', help='List available template names and paths.')
    inspect = commands.add_parser('inspect', help='Read a model without executing it.')
    inspect.add_argument('model', type=Path)
    for name in ('simulate', 'experiment'):
        command = commands.add_parser(name)
        command.add_argument('model', type=Path)
        command.add_argument('--output', type=Path, required=True, help='New output directory.')
        command.add_argument('--seed', type=int)
        if name == 'experiment':
            command.add_argument('--settings', type=Path, help='Sweep settings JSON instead of model.experiments.')
    return parser


def experiment_request(request, settings, seed):
    from gui.experiments import ExperimentRequest, SweepAxis

    def axis(key):
        options = settings.get(key)
        if not isinstance(options, dict) or options.get('target') is None:
            raise ValueError(f'Choose a target for experiments.{key}.')
        target = options['target']
        kind = options.get('kind', 'parameter')
        label = str(target) if kind == 'parameter' else f'Y[{target}]'
        return SweepAxis(kind, target, float(options['minimum']), float(options['maximum']),
                         options.get('count', 10), options.get('scale', 'linear'), label)

    observable = settings.get('observable')
    if observable is None:
        observable = 0
    return ExperimentRequest(request, axis('x'), axis('y'), observable, f'Y[{observable}]',
                             settings.get('metric', 'final_mean'), settings.get('repeats', 1),
                             seed, float(settings.get('tail_percent', 25)) / 100)


def export_simulation(result, path):
    count = result.history[0].shape[1]
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        rows = csv.writer(handle)
        rows.writerow(['frame', 't', 'cell_id', 'x', 'y', *(f'Y[{i}]' for i in range(count))])
        for frame, (values, positions, ids) in enumerate(zip(result.history, result.cell_positions, result.cell_ids, strict=True)):
            for cell_id, point, state in zip(ids, positions, values, strict=True):
                rows.writerow([frame, frame * float(result.params['dT']), cell_id, *point, *state])


def run_job(args, config):
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    from gui.simulation_session import SimulationSession

    output = args.output.resolve()
    if output.exists():
        raise ValueError('Choose a new output directory.')
    settings = {}
    if args.command == 'experiment':
        settings = (json.loads(args.settings.read_text(encoding='utf-8-sig'))
                    if args.settings else config.get('experiments'))
        if not isinstance(settings, dict):
            raise ValueError('Set experiments in the model or supply --settings with a sweep JSON object.')
    seed = args.seed if args.seed is not None else settings.get('seed', 0)
    if type(seed) is not int or not 0 <= seed < 2**32:
        raise ValueError('Seed must be an integer from 0 to 4294967295.')
    np.random.seed(seed)
    grid = model_grid(config)
    initial_positions = grid.copy_cells().tolist()
    request = model_request(config, grid)
    output.parent.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory(prefix='.voronoi-', dir=output.parent) as staging:
        job = Path(staging, 'result')
        job.mkdir()
        if args.command == 'simulate':
            warnings = []
            result = SimulationSession(request).run(warning_callback=warnings.append)
            if any(not np.isfinite(frame).all() for frame in result.history):
                raise ValueError('Simulation contains non-finite values.')
            export_simulation(result, job / 'cells.csv')
            details = {'frames': len(result.history), 'final_cells': len(result.history[-1]), 'warnings': warnings}
        else:
            from matplotlib.figure import Figure
            from gui.experiments import ExperimentRunner
            from gui.experiment_output import draw_heatmap, export_experiment_data

            result = ExperimentRunner(experiment_request(request, settings, seed)).run()
            export_experiment_data(result, job / 'heatmap.csv')
            figure = Figure(figsize=(8, 6), layout='constrained')
            draw_heatmap(figure, result)
            with matplotlib.rc_context({'pdf.fonttype': 42}):
                figure.savefig(job / 'heatmap.pdf')
            details = {'shape': list(result.values.shape),
                       'failed_runs': sum(bool(record['error']) for record in result.records),
                       'warnings': result.warnings}
        write_json(job / 'model.json', config)
        write_json(job / 'initial_grid.json', initial_positions)
        files = sorted(item.name for item in job.iterdir())
        if output.exists():
            raise FileExistsError(f'Output directory already exists: {output}')
        job.rename(output)

    return {'output': str(output), 'seed': seed,
            'files': [str(output / name) for name in files], **details}


def main(argv):
    # Read the report destination before full parsing so syntax errors can be reported.
    report_parser = CliParser(add_help=False)
    report_parser.add_argument('--report', type=Path)
    report = None
    try:
        report = report_parser.parse_known_args(argv)[0].report
        args = cli_parser().parse_args(argv)
        if report:
            inputs = [getattr(args, name, None) for name in ('model', 'settings')]
            if any(path and report.resolve() == path.resolve() for path in inputs):
                report = None
                raise ValueError('Report path must differ from the model and settings files.')
            if args.command in ('simulate', 'experiment') and (
                report.resolve() == args.output.resolve() or args.output.resolve() in report.resolve().parents
            ):
                report = None
                raise ValueError('Keep the report path outside the job output directory.')
        if not cli_enabled():
            raise ValueError('Enable Agent CLI in the application menu before using automation.')
        if args.command == 'list-templates':
            data = {'templates': [{'name': path.stem, 'path': str(path.resolve())} for path in template_files()]}
        else:
            config = read_model(args.model)
            if args.command == 'inspect':
                data = {'model': config}
            else:
                with redirect_stdout(sys.stderr):
                    data = run_job(args, config)
        payload = {'status': 'ok', 'command': args.command, **data}
        code = 0
    except Exception as error:
        payload = {'status': 'error', 'error': str(error)}
        code = 1
    try:
        if report:
            report.parent.mkdir(parents=True, exist_ok=True)
            write_json(report, payload)
        if sys.stdout is not None:
            print(json.dumps(payload, ensure_ascii=True, allow_nan=False))
    except (OSError, ValueError) as error:
        if sys.stderr is not None:
            print(f'Unable to write CLI report: {error}', file=sys.stderr)
        return 1
    return code
