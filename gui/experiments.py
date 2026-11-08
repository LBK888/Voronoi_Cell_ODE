"""Reproducible two-axis sweeps using the application's real simulation engine."""

from dataclasses import asdict, dataclass, replace
from numbers import Real

import numpy as np

from gui.simulation_session import SimulationRequest, SimulationSession
from voronoi_grid import VoronoiGrid


METRICS = {
    'final_mean': 'Final spatial mean',
    'final_std': 'Final spatial standard deviation',
    'final_cv': 'Final spatial CV (std / |mean|)',
    'tail_mean': 'Tail mean of spatial mean',
    'tail_amplitude': 'Tail peak-to-peak of spatial mean',
}


class ExperimentCancelled(Exception):
    pass


@dataclass(frozen=True)
class SweepAxis:
    kind: str  # parameter or initial
    target: object  # parameter name or state-column index
    minimum: float
    maximum: float
    count: int = 10
    scale: str = 'linear'
    label: str = ''

    def values(self):
        if self.kind not in ('parameter', 'initial'):
            raise ValueError('Axis type must be parameter or initial concentration')
        if self.scale not in ('linear', 'log'):
            raise ValueError('Axis scale must be linear or log')
        if not isinstance(self.count, int) or not 2 <= self.count <= 200:
            raise ValueError('Each axis needs 2–200 samples')
        if not np.isfinite([self.minimum, self.maximum]).all() or self.maximum <= self.minimum:
            raise ValueError('Axis maximum must be finite and greater than minimum')
        if self.scale == 'log' and self.minimum <= 0:
            raise ValueError('Logarithmic axes require positive bounds')
        function = np.geomspace if self.scale == 'log' else np.linspace
        return function(self.minimum, self.maximum, self.count)


@dataclass(frozen=True)
class ExperimentRequest:
    simulation: SimulationRequest
    x: SweepAxis
    y: SweepAxis
    observable: int = 0
    observable_label: str = 'Y[0]'
    metric: str = 'final_mean'
    repeats: int = 1
    seed: int = 0
    tail_fraction: float = 0.25


@dataclass
class ExperimentResult:
    request: ExperimentRequest
    x_values: np.ndarray
    y_values: np.ndarray
    values: np.ndarray  # rows = Y, columns = X
    std: np.ndarray  # population standard deviation across replicates
    records: list
    warnings: list

    def metadata(self):
        request = self.request
        return {
            'x': asdict(request.x), 'y': asdict(request.y),
            'observable': request.observable,
            'observable_label': request.observable_label,
            'metric': request.metric, 'repeats': request.repeats,
            'seed': request.seed, 'tail_fraction': request.tail_fraction,
            'T': request.simulation.duration,
            'params': request.simulation.params_source,
            'ode': request.simulation.ode_source,
            'movement': asdict(request.simulation.movement),
            'lifecycle': asdict(request.simulation.lifecycle),
            'parameter_overrides': request.simulation.parameter_overrides,
            'initial_overrides': request.simulation.initial_overrides,
            'initial_cell_positions': request.simulation.grid.copy_cells().tolist(),
            'cell_dist': request.simulation.grid.cell_dist,
            'warnings': self.warnings,
            'time_sampling': 'np.arange(0, T, dT); final frame is before T',
            'aggregation': 'All replicates must succeed; otherwise the pixel is NaN. std uses ddof=0.',
        }


def measure_response(history, observable, metric, tail_fraction):
    """Spatial statistics include all living cells, including boundary cells."""
    for frame in history:
        if not np.isfinite(np.asarray(frame, dtype=float)).all():
            raise ValueError('Simulation contains non-finite values; reduce dT or review the model')
    if not history or not len(history[-1]):
        raise ValueError('No living cells at the final frame')
    last = np.asarray(history[-1], dtype=float)[:, observable]
    if metric == 'final_mean':
        return float(np.mean(last))
    if metric == 'final_std':
        return float(np.std(last))
    if metric == 'final_cv':
        mean = abs(float(np.mean(last)))
        if mean <= 1e-12:
            raise ValueError('CV is undefined because the spatial mean is approximately zero')
        return float(np.std(last) / mean)
    start = max(0, len(history) - max(1, int(np.ceil(len(history) * tail_fraction))))
    tail = history[start:]
    if any(not len(frame) for frame in tail):
        raise ValueError('No living cells in the analysis window')
    means = np.array([np.mean(np.asarray(frame)[:, observable]) for frame in tail])
    if metric == 'tail_mean':
        return float(means.mean())
    if metric == 'tail_amplitude':
        if len(means) < 2:
            raise ValueError('Amplitude requires at least two frames in the tail window')
        return float(np.ptp(means))
    raise ValueError('Unknown response metric')


class ExperimentRunner:
    def __init__(self, request):
        self.request = request
        self.x_values = request.x.values()
        self.y_values = request.y.values()
        if (request.x.kind, request.x.target) == (request.y.kind, request.y.target):
            raise ValueError('X and Y must sweep different targets')
        if request.metric not in METRICS:
            raise ValueError('Choose a valid response metric')
        if not isinstance(request.repeats, int) or not 1 <= request.repeats <= 100:
            raise ValueError('Replicates must be between 1 and 100')
        if not isinstance(request.seed, int) or not 0 <= request.seed <= 2**32 - 1:
            raise ValueError('Seed must be an integer between 0 and 2**32 - 1')
        if not 0 < request.tail_fraction <= 1:
            raise ValueError('Tail fraction must be in (0, 1]')

    def run(self, progress_callback=None, check_cancel=None):
        request = self.request
        random_state = np.random.get_state()
        try:
            np.random.seed(request.seed)
            baseline = SimulationSession(request.simulation)
            n_var = baseline.init_y.shape[1]
            if not isinstance(request.observable, int) or not 0 <= request.observable < n_var:
                raise ValueError('Observable is outside the state vector')
            for axis in (request.x, request.y):
                if axis.kind == 'initial':
                    if not isinstance(axis.target, int) or not 0 <= axis.target < n_var:
                        raise ValueError('Initial concentration target is outside the state vector')
                else:
                    value = baseline.params.get(axis.target)
                    if (not isinstance(value, Real) or isinstance(value, bool)
                            or axis.target == 'n_var' or not np.isfinite(value)):
                        raise ValueError(f'{axis.target} must be an existing finite scalar model parameter')
                if axis.kind == 'parameter' and axis.target == 'dT':
                    if axis.minimum <= 0 or axis.maximum >= request.simulation.duration:
                        raise ValueError('Swept dT must be positive and smaller than T')
            if not np.isfinite(request.simulation.duration) or baseline.params['dT'] >= request.simulation.duration:
                raise ValueError('T must be finite and greater than dT for an experiment')
            return self._sweep(progress_callback, check_cancel)
        finally:
            np.random.set_state(random_state)

    def _sweep(self, progress_callback, check_cancel):
        request = self.request
        shape = (len(self.y_values), len(self.x_values))
        values, deviations = np.full(shape, np.nan), np.full(shape, np.nan)
        records, warnings = [], []
        total = shape[0] * shape[1] * request.repeats
        completed = 0
        last_progress = -1

        def report(percent, message):
            nonlocal last_progress
            if progress_callback and percent != last_progress:
                last_progress = percent
                progress_callback(percent, message)

        def ensure_running():
            if check_cancel and check_cancel():
                raise ExperimentCancelled()

        def warn(message):
            if message not in warnings and len(warnings) < 20:
                warnings.append(message)

        for row, y in enumerate(self.y_values):
            for column, x in enumerate(self.x_values):
                responses = []
                for repeat in range(request.repeats):
                    ensure_running()
                    # Common random numbers: the same replicate starts with the
                    # same seed at every pixel, making parameter effects comparable.
                    seed = (request.seed + repeat) % 2**32
                    np.random.seed(seed)
                    overrides = dict(request.simulation.parameter_overrides)
                    initial = dict(request.simulation.initial_overrides)
                    for axis, value in ((request.x, x), (request.y, y)):
                        (overrides if axis.kind == 'parameter' else initial)[axis.target] = float(value)
                    response, error = np.nan, ''
                    try:
                        simulation = replace(
                            request.simulation,
                            grid=VoronoiGrid.from_cells(request.simulation.grid.copy_cells(),
                                                       cell_dist=request.simulation.grid.cell_dist),
                            parameter_overrides=overrides, initial_overrides=initial,
                        )
                        session = SimulationSession(simulation)

                        def progress(step, steps):
                            ensure_running()
                            percent = int(100 * (completed + step / max(steps, 1)) / total)
                            report(percent, f'Run {completed + 1}/{total}: X={x:.4g}, Y={y:.4g}')

                        result = session.run(progress_callback=progress, check_cancel=check_cancel,
                                             warning_callback=warn)
                        ensure_running()
                        response = measure_response(result.history, request.observable,
                                                    request.metric, request.tail_fraction)
                    except ExperimentCancelled:
                        raise
                    except Exception as exception:
                        error = str(exception)
                    records.append(dict(x=float(x), y=float(y), replicate=repeat + 1,
                                        seed=seed, response=response, error=error))
                    responses.append(response)
                    completed += 1
                    report(int(100 * completed / total), f'Completed {completed}/{total} runs')
                if np.isfinite(responses).all():
                    values[row, column] = np.mean(responses)
                    deviations[row, column] = np.std(responses)
        ensure_running()
        return ExperimentResult(request, self.x_values, self.y_values,
                                values, deviations, records, warnings)
