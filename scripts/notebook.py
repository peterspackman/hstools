from .decompose import shift_to_origin
import pandas as pd
import numpy as np
import sbf

import plotly.offline as py
import plotly.graph_objs as go
py.init_notebook_mode()

def read_hs(filename, with_colors='d_norm'):
    f = sbf.File(filename)
    f.read()
    vertices = f['vertices'].data.transpose()
    vertices = shift_to_origin(vertices)
    faces = f['faces'].data.transpose() -1
    colors = with_colors, f[with_colors].data
    return vertices, faces, colors   

def create_figure(surface, colorscale='Viridis'):
    vertices, faces, colors = surface
    data = go.Data([
        go.Mesh3d(
            x = vertices[:, 0],
            y = vertices[:, 1],
            z = vertices[:, 2],
            colorbar = go.ColorBar(title=colors[0]),
            colorscale = colorscale,
            reversescale=True,
            intensity = colors[1],
            i = faces[:, 0],
            j = faces[:, 1],
            k = faces[:, 2],
            name='y',
            showscale=True,
        )
    ])
    axis_kws = dict(
        title='',
        showaxeslabels=False,
        showbackground=False,
        zeroline=False,
        ticks='',
        showticklabels=False,
        showgrid=False,
    )
    layout = go.Layout(
        scene = dict(
            xaxis = axis_kws,
            yaxis = axis_kws,
            zaxis = axis_kws,
        ),
    )
    return go.Figure(data=data, layout=layout)

def show_surface(fig):
    return py.iplot(fig)

def plot_hs(filename, with_colors='d_norm', colorscale='Viridis'):
    surface = read_hs(filename, with_colors=with_colors)
    fig = create_figure(surface, colorscale=colorscale)
    return show_surface(fig)
