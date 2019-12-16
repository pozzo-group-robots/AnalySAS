import numpy as np
import plotly.graph_objects as go
import ipywidgets as widgets


def plot_sas(sas_data):

    """
    Returns a widget object to plot one or more SasData class objects.

    Parameters:
    -----------

    sas_data: SasData or iterable of SasData
        Scattering data provided as a single SasData object or iterable of SasData objects.

    """
    try:
        sas_data = list(sas_data)
    except:
        sas_data = [sas_data]
    else:
        pass

    keys = [x for x in range(0,len(sas_data))]
    data_dict = dict(zip(keys, sas_data))
    names = [x.label for x in sas_data]
    name_dict = dict(zip(keys, names))

    # Widget Components

    select = widgets.SelectMultiple(
        options=[(v, k) for k, v in name_dict.items()],
        description='Samples (select multiple):',
        style = {'description_width': 'initial'},
        layout=widgets.Layout(flex='1 1 0%', width='auto', height='120px'),
        disabled=False)

    update_button = widgets.Button(
        description='Update Plot',
        disabled=False,
        layout=widgets.Layout(flex='1 1 0%', width='auto'),
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Update Plot')

    add_button = widgets.Button(
        description='Add to Plot',
        disabled=False,
        layout=widgets.Layout(flex='1 1 0%', width='auto'),
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Add to Plot')

    clear_button = widgets.Button(
        description='Clear Plot',
        disabled=False,
        layout=widgets.Layout(flex='1 1 0%', width='auto'),
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Clear Plot')

    buttons = widgets.VBox([update_button, add_button, clear_button])
    selection = widgets.HBox([select, buttons])

    fig = go.FigureWidget()
    fig = fig.update_layout(xaxis_type="log", yaxis_type="log", height=500, showlegend=True, margin=dict(t=30))

    fig = fig.update_xaxes(exponentformat='e', showexponent='all')
    fig = fig.update_yaxes(exponentformat='e', showexponent='all')

    widget = widgets.VBox([selection,fig])
    widget

    # Interactive Widget Functions

    def update_plot(b):
        fig.data = []
        for value in select.value:
            data = data_dict[value]
            fig.add_scatter(x = data.q, y = data.Iq, error_y = dict(array=data.dIq), mode='markers', name=data.label)

    def add_to_plot(b):
        current_data = []
        for trace in fig.data:
            current_data.append(trace['name'])
        for value in select.value:
            data = data_dict[value]
            if data.label not in current_data:
                fig.add_scatter(x = data.q, y = data.Iq, error_y = dict(array=data.dIq), mode='markers', name=data.label)

    def clear_plot(b):
        fig.data = []

    # Define Widget Responses

    update_button.on_click(update_plot)
    add_button.on_click(add_to_plot)
    clear_button.on_click(clear_plot)

    return(widget)
