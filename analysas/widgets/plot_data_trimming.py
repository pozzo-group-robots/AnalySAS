import numpy as np
import plotly.graph_objects as go
import ipywidgets as widgets


def plot_data_trimming(sas_data):

    """
    Returns a widget object to interactively trim scattering data stored in a SasData object.

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
    button_prev = widgets.Button(
        description='Previous',
        layout=widgets.Layout(flex='0.5 1 0%', width='auto'),
        button_style='',
        tooltip='Previous')

    dropdown = widgets.Dropdown(
        options=[(v, k) for k, v in name_dict.items()],
        description='Sample:',
        layout=widgets.Layout(flex='2 1 0%', width='auto'),
        disabled=False)

    button_next = widgets.Button(
        description='Next',
        layout=widgets.Layout(flex='0.5 1 0%', width='auto'),
        button_style='',
        tooltip='Next')

    header = widgets.HBox([button_prev, dropdown, button_next])


    q_min = widgets.Text(
        value='',
        placeholder='',
        description='Min Q:',
        layout=widgets.Layout(flex='0.25 1 0%', width='auto'),
        disabled=False)

    q_max = widgets.Text(
        value='',
        placeholder='',
        description='Max Q:',
        layout=widgets.Layout(flex='0.25 1 0%', width='auto'),
        disabled=False)

    ind_vals = widgets.Text(
        value='',
        placeholder='',
        description='Trim Points:',
        layout=widgets.Layout(flex='1 1 0%', width='auto'),
        disabled=False)

    trim_info = widgets.HBox([q_min, q_max, ind_vals])

    button_reset = widgets.Button(
        description='Reset Trim',
        button_style='',
        layout=widgets.Layout(flex='0.5 1 0%', width='auto'),
        tooltip='Reset Trim')

    button_trim = widgets.Button(
        description='Trim Data',
        button_style='',
        layout=widgets.Layout(flex='0.5 1 0%', width='auto'),
        tooltip='Trim Data')

    trim_buttons = widgets.HBox([button_reset, button_trim])

    data = data_dict[dropdown.value]
    data_plot = go.Scatter(x = data.q, y = data.Iq, error_y = dict(array=data.dIq), mode='markers')

    fig = go.FigureWidget(data=[data_plot])
    fig = fig.update_layout(xaxis_type="log", yaxis_type="log", title=data.label, height=500)

    fig = fig.update_xaxes(exponentformat='e', showexponent='all')
    fig = fig.update_yaxes(exponentformat='e', showexponent='all')

    widget = widgets.VBox([header, trim_info, trim_buttons, fig])

    # Interactive Widget Functions

    def on_prev_click(b):
        current_val = dropdown.value
        values = np.array([k for (v, k) in list(dropdown.options)])
        current_ind = np.where(values == current_val)
        if current_ind[0] == 0:
            pass
        else:
            dropdown.value = values[current_ind[0] - 1]


    def on_next_click(b):
        current_val = dropdown.value
        values = np.array([k for (v, k) in list(dropdown.options)])
        current_ind = np.where(values == current_val)
        if current_ind[0] == len(values) - 1:
            pass
        else:
            dropdown.value = values[current_ind[0] + 1]

    def update_plot(d):
        data = data_dict[dropdown.value]

        with fig.batch_update():
                fig.data[0].x, fig.data[0].y, fig.data[0].error_y = data.q, data.Iq, dict(array=data.dIq)

        fig.update_layout(title=data.label)

    def update_trim_boxes(d):
        q_min.value = ''
        q_max.value = ''
        ind_vals.value = ''

    def on_reset_click(b):
        data = data_dict[dropdown.value]
        data.reset_trim()

    def on_trim_click(b):
        data = data_dict[dropdown.value]
        data.reset_trim()

        # formatting q_min
        if len(q_min.value) > 0:
            set_min = float(q_min.value)
        else:
            set_min = np.min(data.q)/5

        # formatting q_max
        if len(q_max.value) > 0:
            set_max = float(q_max.value)
        else:
            set_max = np.max(data.q)*5

        data.trim_q_range(set_min, set_max)

        if len(ind_vals.value) > 0:
            set_vals = ind_vals.value
            set_vals = set_vals.replace(' ','')
            set_vals = set_vals.split(',')

            ints = [int(x) for x in set_vals if (x[0] is not '(' and x[-1] is not ')')]
            lows = [int(x[1:]) for x in set_vals if x[0] is '(']
            highs = [int(x[:-1]) for x in set_vals if x[-1] is ')']
            tuples = [(x,y) for x, y in zip(lows, highs)]
            trim_points = ints + tuples

            data.trim_points(trim_points)


     # Define Widget Responses

    button_prev.on_click(on_prev_click)
    button_next.on_click(on_next_click)
    dropdown.observe(update_plot, names='value')
    dropdown.observe(update_trim_boxes, names='value')
    button_reset.on_click(on_reset_click)
    button_reset.on_click(update_trim_boxes)
    button_reset.on_click(update_plot)
    button_trim.on_click(on_trim_click)
    button_trim.on_click(update_plot)

    return(widget)
