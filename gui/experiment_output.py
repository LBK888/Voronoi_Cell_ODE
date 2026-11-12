"""Write and render sweep results for both the GUI and the CLI."""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import colormaps

from gui.experiments import METRICS
from gui.model_files import write_json


def _edges(values, scale):
    transformed = np.log(values) if scale == 'log' else values
    edges = np.r_[transformed[0] - (transformed[1] - transformed[0]) / 2,
                  (transformed[:-1] + transformed[1:]) / 2,
                  transformed[-1] + (transformed[-1] - transformed[-2]) / 2]
    return np.exp(edges) if scale == 'log' else edges


def draw_heatmap(figure, result):
    figure.clear()
    axis = figure.add_subplot(111)
    request = result.request
    cmap = colormaps['viridis'].with_extremes(bad='#C9D1D9')
    image = axis.pcolormesh(_edges(result.x_values, request.x.scale),
                           _edges(result.y_values, request.y.scale),
                           np.ma.masked_invalid(result.values), cmap=cmap, shading='flat')
    axis.set_xscale(request.x.scale)
    axis.set_yscale(request.y.scale)
    axis.set_xlabel(request.x.label + (' (initial)' if request.x.kind == 'initial' else ''))
    axis.set_ylabel(request.y.label + (' (initial)' if request.y.kind == 'initial' else ''))
    metric = METRICS[request.metric]
    axis.set_title(f'{request.observable_label}\n{metric}')
    figure.colorbar(image, ax=axis, label=metric)


def export_experiment_data(result, path):
    path = Path(path)
    table = pd.DataFrame(result.records)
    stats = {
        (x, y): (result.values[row, col], result.std[row, col])
        for row, y in enumerate(result.y_values) for col, x in enumerate(result.x_values)
    }
    table[['pixel_mean', 'pixel_std']] = [stats[(r.x, r.y)] for r in table.itertuples()]
    table.to_csv(path, index=False, encoding='utf-8-sig')
    write_json(path.with_name(path.stem + '_settings.json'), result.metadata())
