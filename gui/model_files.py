"""Model JSON and external template locations shared by the GUI and CLI."""

import ast
import json
import math
import os
from pathlib import Path
import sys
import tempfile

from gui.app_paths import data_path, resource_path


GRID_MODES = ('honeycomb', 'random', 'regular', 'import')
LIFECYCLE_MODES = ('area', 'random')
MOVEMENT_OPTIONS = (
    ('random', 'random_strength', '0.05'),
    ('away', 'away_strength', '0.1'),
    ('ce', 'convergent_strength', '0.02'),
    ('repulsion', 'repulsion_strength', '0.08'),
    ('vertex_model', 'vertex_strength', '0.05'),
)


def user_templates():
    return data_path('templates')


def template_files():
    external = Path(sys.executable).parent / 'templates' if getattr(sys, 'frozen', False) else resource_path('templates')
    files = {}
    for directory in (external, user_templates()):
        if directory.is_dir():
            for path in sorted(directory.iterdir()):
                if path.is_file() and path.suffix.lower() == '.json':
                    files[path.name.casefold()] = path
    return sorted(files.values(), key=lambda path: path.name.casefold())


def write_json(path, value):
    """Replace a JSON file only after serialization and writing succeed."""
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                         delete=False, suffix='.tmp') as handle:
            temporary = Path(handle.name)
            json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write('\n')
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def read_model(path):
    path = Path(path)
    if path.suffix.lower() != '.json':
        raise ValueError('Choose a model JSON file.')
    config = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(config, dict) or not {'ode', 'params'}.intersection(config):
        raise ValueError('A model JSON must contain ode or params.')
    return normalize_model(config)


def normalize_model(config):
    """Validate widget inputs and Python syntax without executing model code."""
    from gui.default_color_func import get_default_color_func
    from gui.default_ode import get_default_ode
    from gui.default_params import get_default_params

    config = dict(config)
    for key, default in (
        ('ode', get_default_ode), ('params', get_default_params),
        ('color_func', get_default_color_func),
    ):
        value = config.get(key, '')
        if not isinstance(value, str):
            raise ValueError(f'{key} must be text.')
        value = value.strip() or default()
        ast.parse(value, filename=key)
        config[key] = value

    integer_fields = {
        'grid_shape_x': (20, 1, 200), 'grid_shape_y': (20, 1, 200),
        'grid_mode': (0, 0, 3),
        'move_division_n': (1, 1, 10), 'move_apoptosis_n': (1, 1, 10),
        'move_division_method': (0, 0, 1), 'move_apoptosis_method': (0, 0, 1),
    }
    for key, (default, low, high) in integer_fields.items():
        value = config.get(key, default)
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f'{key} must be an integer from {low} to {high}.')
        config[key] = value

    for name in ('division', 'apoptosis', *(option[0] for option in MOVEMENT_OPTIONS)):
        key = f'move_{name}'
        value = config.get(key, False)
        if type(value) is not bool:
            raise ValueError(f'{key} must be true or false.')
        config[key] = value

    numeric_fields = {'cell_dist': '1.0', 'pos_rand': '0.25', 'T': '30.0'}
    numeric_fields.update({f'move_{name}_strength': default for name, _, default in MOVEMENT_OPTIONS})
    for key, default in numeric_fields.items():
        value = config.get(key, default)
        try:
            number = float(value)
        except (TypeError, ValueError):
            raise ValueError(f'{key} must be a finite number.') from None
        if isinstance(value, bool) or not math.isfinite(number):
            raise ValueError(f'{key} must be a finite number.')
        if key in ('cell_dist', 'T') and number <= 0:
            raise ValueError(f'{key} must be greater than zero.')
        config[key] = str(value)
    for key in ('replicate', 'repeats'):
        config[key] = str(config.get(key, '5'))

    if 'cell_positions' in config:
        points = config['cell_positions']
        if not isinstance(points, list) or len(points) < 3 or any(
            not isinstance(point, list) or len(point) != 2
            or any(type(value) not in (int, float) or not math.isfinite(value) for value in point)
            for point in points
        ):
            raise ValueError('cell_positions must contain at least three finite [x, y] pairs.')
    return config


def model_grid(config):
    from voronoi_grid import VoronoiGrid

    mode = GRID_MODES[config['grid_mode']]
    distance = float(config['cell_dist'])
    if mode == 'import':
        if 'cell_positions' not in config:
            raise ValueError('Import the cell coordinates in the GUI and save a template before running this grid from the CLI.')
        return VoronoiGrid.from_cells(config['cell_positions'], cell_dist=distance)
    return VoronoiGrid(grid_shape=(config['grid_shape_x'], config['grid_shape_y']),
                       cell_dist=distance, pos_rand=float(config['pos_rand']), mode=mode)


def model_request(config, grid):
    from gui.simulation_session import CellLifecyclePlan, SimulationRequest
    from voronoi_grid import MovementPlan

    movement = MovementPlan(**{
        field: float(config[f'move_{name}_strength']) if config[f'move_{name}'] else 0.0
        for name, field, _ in MOVEMENT_OPTIONS
    })
    lifecycle = CellLifecyclePlan(
        division_count=config['move_division_n'] if config['move_division'] else 0,
        division_mode=LIFECYCLE_MODES[config['move_division_method']],
        apoptosis_count=config['move_apoptosis_n'] if config['move_apoptosis'] else 0,
        apoptosis_mode=LIFECYCLE_MODES[config['move_apoptosis_method']],
    )
    return SimulationRequest(grid, float(config['T']), config['params'], config['ode'], movement, lifecycle)
