"""Prepare and run one simulation behind a small, testable interface."""

from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import numpy as np

from biophysics_model import BiophysicsModel
from gui.sim_utils import diffusion_weighted_mean, get_voronoi_neighbors_wo_outer
from voronoi_grid import MovementPlan, VoronoiGrid


class SimulationConfigurationError(ValueError):
    pass


@dataclass(frozen=True)
class CellLifecyclePlan:
    division_count: int = 0
    division_mode: str = 'area'
    apoptosis_count: int = 0
    apoptosis_mode: str = 'area'


@dataclass(frozen=True)
class SimulationRequest:
    grid: VoronoiGrid
    duration: float
    params_source: str
    ode_source: str
    movement: MovementPlan = field(default_factory=MovementPlan)
    lifecycle: CellLifecyclePlan = field(default_factory=CellLifecyclePlan)


@dataclass(frozen=True)
class SimulationResult:
    history: list
    cell_positions: list
    grid: VoronoiGrid
    params: dict
    cell_ids: list = field(default_factory=list)


class SimulationSession:
    """Deep module that owns simulation preparation, scheduling, and execution."""

    def __init__(self, request):
        if not isinstance(request, SimulationRequest):
            raise TypeError('request must be a SimulationRequest')
        if not isinstance(request.grid, VoronoiGrid):
            raise TypeError('request.grid must be a VoronoiGrid')
        if request.duration <= 0:
            raise SimulationConfigurationError('duration must be positive')

        self.request = request
        self.params = self._parse_params(request.params_source)
        self.ode_func, self.mechanics_func = self._compile_model(request.ode_source)
        self.init_y = self._build_initial_state(len(request.grid.cells))
        self.model = BiophysicsModel(
            cell_count=len(request.grid.cells),
            params=self.params,
            ode_func=self.ode_func,
            mechanics_func=self.mechanics_func,
            init_Y=self.init_y,
            vor_grid=request.grid,
            movement_plan=request.movement,
        )

    def run(self, progress_callback=None, check_cancel=None,
            warning_callback=None):
        lifecycle = self.request.lifecycle
        scheduled_steps = self._scheduled_steps()
        history, cell_positions = self.model.simulate(
            T=self.request.duration,
            proliferation_steps=scheduled_steps if lifecycle.division_count else [],
            apoptosis_steps=scheduled_steps if lifecycle.apoptosis_count else [],
            proliferation_n=lifecycle.division_count,
            apoptosis_n=lifecycle.apoptosis_count,
            proliferation_mode=lifecycle.division_mode,
            apoptosis_mode=lifecycle.apoptosis_mode,
            progress_callback=progress_callback,
            check_cancel=check_cancel,
            warning_callback=warning_callback,
        )
        return SimulationResult(
            history=list(history),
            cell_positions=list(cell_positions),
            grid=self.model.vor_grid,
            params=self.params,
            cell_ids=list(self.model.cell_ids_history),
        )

    def _parse_params(self, source):
        expression = source.strip()
        if expression.startswith('params ='):
            expression = expression[len('params ='):].strip()
        try:
            params = eval(
                expression,
                {'__builtins__': {}, 'dict': dict},
                {},
            )
        except Exception as error:
            raise SimulationConfigurationError(
                f'Unable to parse params: {error}'
            ) from error
        if not isinstance(params, dict):
            raise SimulationConfigurationError('params must evaluate to a dict')
        try:
            d_t = float(params['dT'])
        except (KeyError, TypeError, ValueError) as error:
            raise SimulationConfigurationError(
                "params['dT'] must be a positive number"
            ) from error
        if d_t <= 0:
            raise SimulationConfigurationError("params['dT'] must be positive")
        return params

    @staticmethod
    def _compile_model(source):
        globals_dict = {
            'np': np,
            'plt': plt,
            'get_voronoi_neighbors_wo_outer': get_voronoi_neighbors_wo_outer,
            'diffusion_weighted_mean': diffusion_weighted_mean,
        }
        locals_dict = {}
        try:
            exec(source, globals_dict, locals_dict)
        except Exception as error:
            raise SimulationConfigurationError(
                f'Unable to compile ODE: {error}'
            ) from error

        ode_func = locals_dict.get('ode')
        mechanics_func = locals_dict.get('mechanics')
        if not callable(ode_func):
            raise SimulationConfigurationError(
                'ODE source must define a callable ode function'
            )
        return ode_func, mechanics_func

    def _build_initial_state(self, cell_count):
        if 'init_Y' in self.params:
            initial = np.asarray(self.params['init_Y'], dtype=float)
            if initial.ndim == 1:
                initial = np.tile(initial, (cell_count, 1))
            elif initial.ndim == 2 and initial.shape[0] == 1:
                initial = np.repeat(initial, cell_count, axis=0)
            if initial.ndim != 2 or initial.shape[0] != cell_count:
                raise SimulationConfigurationError(
                    'init_Y must be a state vector or have one row per cell'
                )
            self.params['n_var'] = initial.shape[1]
            return initial.copy()

        try:
            n_variables = int(self.params.get('n_var', 4))
            sigma = float(self.params.get('sigma', 0.0) or 0.0)
        except (TypeError, ValueError) as error:
            raise SimulationConfigurationError(
                'n_var and sigma must be numeric'
            ) from error
        if n_variables <= 0:
            raise SimulationConfigurationError('n_var must be positive')
        return np.random.randn(cell_count, n_variables) * sigma

    def _scheduled_steps(self):
        d_t = float(self.params['dT'])
        total_steps = int(self.request.duration / d_t)
        step_interval = max(1, int(1 / d_t))
        return list(range(0, total_steps, step_interval))
