"""Validate a trusted model JSON against the actual GUI simulation contract."""

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys


def validate(config_path, repo, duration=None, seed=0, halving=False, rtol=.05):
    sys.path.insert(0, str(repo))
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba_array
    import numpy as np
    from gui.simulation_session import SimulationRequest, SimulationSession
    from voronoi_grid import VoronoiGrid

    config = json.loads(config_path.read_text(encoding='utf-8'))
    grid = VoronoiGrid(grid_shape=(5, 5), pos_rand=0)
    run_time = float(duration if duration is not None else min(float(config.get('T', 1)), 1))
    request = SimulationRequest(grid, run_time, config['params'], config['ode'])
    np.random.seed(seed)
    session = SimulationSession(request)
    if run_time <= session.params['dT']:
        raise ValueError('Validation duration must exceed dT')
    state = session.init_y.copy()
    snapshot = state.copy()
    derivative = np.asarray(session.ode_func(state, 0., session.params, grid.cells.copy()))
    if derivative.shape != state.shape or not np.isfinite(derivative).all():
        raise ValueError('ode must return finite derivatives of the same shape as Y')
    if not np.array_equal(state, snapshot):
        raise ValueError('ode modifies its input state')
    labels = session.params.get('labels', [])
    if len(labels) != state.shape[1]:
        raise ValueError('Provide exactly one label per state column')
    if config.get('color_func'):
        namespace = {'np': np, 'plt': plt}
        exec(config['color_func'], namespace)
        import inspect
        function = namespace['color_func']
        for mode in ('polygon', 'center', 'off_center', 'off_center2', 'membrane'):
            kwargs = {'mode': mode} if 'mode' in inspect.signature(function).parameters else {}
            if 'global_bounds' in inspect.signature(function).parameters:
                kwargs['global_bounds'] = [(float(v.min()), float(v.max())) for v in state.T]
            colors = function(state, **kwargs)
            if colors is not None and to_rgba_array(colors).shape != (len(state), 4):
                raise ValueError(f'Invalid color array for mode {mode}')
    np.random.seed(seed)
    session = SimulationSession(request)
    result = session.run()
    if any(not np.isfinite(frame).all() for frame in result.history):
        raise ValueError('Non-finite simulation result')
    report = dict(status='passed', cells=len(state), variables=state.shape[1],
                  frames=len(result.history),
                  final_time=(len(result.history) - 1) * session.params['dT'],
                  minimum=float(min(frame.min() for frame in result.history)))
    if halving:
        if session.mechanics_func is not None:
            raise ValueError('Halving comparison requires a fixed-cell deterministic model')
        np.random.seed(seed)
        refined = SimulationSession(replace(
            request, grid=VoronoiGrid(grid_shape=(5, 5), pos_rand=0),
            parameter_overrides={'dT': session.params['dT'] / 2},
        )).run()
        fine = refined.history[2 * (len(result.history) - 1)]
        coarse = result.history[-1]
        difference = float(np.linalg.norm(coarse - fine) / max(np.linalg.norm(fine), 1e-12))
        report['halving_relative_error_at_common_time'] = difference
        if difference > rtol:
            raise ValueError(f'Halving comparison error {difference:.4g} exceeds {rtol}; refine dT')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--duration', type=float)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--halving-check', action='store_true', help='For deterministic, fixed-cell models only')
    parser.add_argument('--rtol', type=float, default=.05)
    args = parser.parse_args()
    try:
        print(json.dumps(validate(args.config, args.repo, args.duration, args.seed,
                                  args.halving_check, args.rtol), indent=2))
    except Exception as error:
        print(f'Model validation failed: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
