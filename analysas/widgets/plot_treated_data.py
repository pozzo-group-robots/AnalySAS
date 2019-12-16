import numpy as np
import plotly.graph_objects as go
import ipywidgets as widgets


def plot_treated_data(sas_data):

    """
    Returns a widget object to inspect the raw and treated data of one or more SasData class objects.

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

    data = data_dict[dropdown.value]
    q_raw, Iq_raw, dIq_raw, dq_raw = data.get_raw_data()
    incoh_bkgd = data.get_incoh_bkgd()
    raw_data = go.Scatter(x = q_raw, y = Iq_raw, error_y=dict(array=dIq_raw), mode='markers', name='Raw Scattering Data')
    inc_back = go.Scatter(x = [np.min(q_raw), np.max(q_raw)], y = [incoh_bkgd, incoh_bkgd], mode='lines', name='Incoherent Background')
    sub_data = go.Scatter(x = data.q, y = data.Iq, error_y = dict(array=data.dIq), mode='markers', name='Treated Scattering Data')

    fig = go.FigureWidget(data=[raw_data, sub_data, inc_back])
    fig = fig.update_layout(xaxis_type="log", yaxis_type="log", title=data.label, height=500)

    fig = fig.update_xaxes(exponentformat='e', showexponent='all')
    fig = fig.update_yaxes(exponentformat='e', showexponent='all')

    widget = widgets.VBox([header, fig])

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
        q_raw, Iq_raw, dIq_raw, dq_raw = data.get_raw_data()
        incoh_bkgd = data.get_incoh_bkgd()

        with fig.batch_update():
                fig.data[0].x, fig.data[0].y, fig.data[0].error_y = q_raw, Iq_raw, dict(array=dIq_raw)
                fig.data[1].x, fig.data[1].y, fig.data[1].error_y = data.q, data.Iq, dict(array=data.dIq)
                fig.data[2].x, fig.data[2].y =[np.min(q_raw),np.max(q_raw)], [incoh_bkgd, incoh_bkgd]

        fig.update_layout(title=data.label)

    # Define Widget Responses

    button_prev.on_click(on_prev_click)
    button_next.on_click(on_next_click)
    dropdown.observe(update_plot, names='value')

    return(widget)
