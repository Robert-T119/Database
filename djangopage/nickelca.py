from djangopage.reactions import *
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import least_squares
from django_plotly_dash import DjangoDash

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app = DjangoDash('nickelca', add_bootstrap_links=True)
# , add_bootstrap_links=True
# add_bootstrap_links=True


app.layout = html.Div([
    html.Div([
        html.H2(u'Nickel-Citrate-Ammonia-H\u2082O system')],
        style={
            'text-align':'center',
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '-10px 0 -10px 0',
            'margin-bottom': '2px',
            'width': '100%',
        }
        ),

    html.Div([
        html.H6(u"Total nickel concentration (kmolm\u207B\u00B3):"),
        dcc.Slider(
            id='nickel_slider',
            min=0.1,
            max=3.0,
            value=1.0,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3]},
                ),
            ],
        style={
            'padding': '10px 15px 10px',
            'border': 'thin lightgrey solid',
            'margin-bottom': '3px',
            'backgroundColor': 'rgb(250, 250, 250)',
            'width': '100%',
            }
        ),

    html.Div([
            html.H6(u"Total citrate concentration (kmolm\u207B\u00B3):"),
            dcc.Slider(
                id='citrate_slider',
                min=0.1,
                max=3.0,
                value=1.1,
                step=0,
                marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3,
                    0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                    1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3]},
                    ),
                ],
            style={
                'padding': '10px 15px 10px',
                'border': 'thin lightgrey solid',
                'margin-bottom': '3px',
                'backgroundColor': 'rgb(250, 250, 250)',
                'width': '100%',
                }
            ),

    html.Div([
        html.H6(u'Total ammonia concentration (kmolm\u207B\u00B3):'),
        dcc.Slider(
            id='ammonia_dropdown',
            min=0.1,
            max=3.0,
            value=1.2,
            step=0,
            marks={i: str(i) for i in [0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3]},
        ),
    ],
        style={
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250,250,250)',
            'padding': '10px 15px 10px',
            'margin-bottom': '3px',
            'width': '100%',
            }
    ),

    html.Div([
        html.Div([
            dcc.Graph(
                # animate=True,
                id='speciation plot',
                )
        ],
        style={
            'width': '48%',
            'margin-right': '1%',
        }
        ),
        html.Div([
            dcc.Graph(
                id='Potential-pH curve',
            ),
        ],
        style={
            'width': '48%',
            'margin-left': '1%',
        }
        )
    ],
    className='row',
        style={
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '0 0px 0 15px',
            "margin-left": "0px",
            'margin-right': '0px',
            'margin-bottom': '3px'
            }),


], style={
        'margin-top': '40px',
        # "margin-left": "10px",
        # 'margin-right': '10px',
        'margin-bottom': '5px'
})

@app.callback(
    Output('speciation plot', 'figure'),
    [Input('nickel_slider', 'value'),
    Input('citrate_slider', 'value'),
    Input('ammonia_dropdown', 'value')])

def speciation_graph(ni_total, citrate_total,ammonia_total):
    k1 = 9.197 * (10 ** 11)
    k2 = 9.905 * (10 ** 2)
    k3 = 5.595 * (10 ** 4)
    k4 = 2.481 * (10 ** 6)
    k5 = 2.512 * (10 ** 5)
    k6 = 1.995 * (10 ** 3)
    k7 = 5.623 * (10 ** 1)
    k8 = 10 ** (-9.25)
    k9 = 10 ** (-18.26)
    k10 = 10 ** (-35.91)
    pH_x = np.linspace(0, 14, 71)
    T_=298
    def concs(citrate_total, ni_total, ammonia_total, pH_x):
        h = 10 ** (-pH_x)
        def equations(p):
            cit3 = p[0]
            nio2 = p[1]
            nh4 = p[2]
            Hcit = h * cit3 * k4
            H2cit = h * Hcit * k3
            H3cit = h * H2cit * k2
            ni2pfree = (nio2 * k1 * (h ** 2)) / (1 + ((nio2 * k1 * (h ** 2)) / ni_total))
            NiH2cit = k7 * ni2pfree * H2cit
            NiHcit = k6 * ni2pfree * Hcit
            Nicit = k5 * ni2pfree * cit3
            nh3 = k8 * nh4 / h
            nin4p2 = k9 * nio2 * (nh4 ** 4) / (h ** 2)
            nin6p2 = k10 * nio2 * (nh4 ** 6) / (h ** 4)
            return (citrate_total - Hcit - H2cit - H3cit - Nicit - NiHcit - NiH2cit - cit3,
                    ni_total - Nicit - NiHcit - NiH2cit - ni2pfree - nin4p2 - nin6p2 - nio2,
                    ammonia_total - nh3 - 4 * nin4p2 - 6 * nin6p2 - nh4)
        res = least_squares(equations, (0.1, 0.1, 0.1), bounds=((0, 0, 0), (3, 3, 3)), method='dogbox', xtol=1e-12)
        # res = least_squares(equations,(0.1,0.1,0.1),method='lm')
        # res = root(equations,[0.1,0.1,0.1],method='lm')
        cit3 = res.x[0]
        nio2 = res.x[1]
        nh4 = res.x[2]
        ni2pfree = (nio2 * k1 * (h ** 2)) / (1 + ((nio2 * k1 * (h ** 2)) / ni_total))
        Hcit = h * cit3 * k4
        H2cit = h * Hcit * k3
        H3cit = h * H2cit * k2
        NiH2cit = k7 * ni2pfree * H2cit
        NiHcit = k6 * ni2pfree * Hcit
        Nicit = k5 * ni2pfree * cit3
        nh3 = k8 * nh4 / h
        nin4p2 = k9 * nio2 * (nh4 ** 4) / (h ** 2)
        nin6p2 = k10 * nio2 * (nh4 ** 6) / (h ** 4)
        return [cit3, nio2, nh4, ni2pfree, Hcit, H2cit, H3cit, NiH2cit, NiHcit, Nicit, nh3, nin4p2, nin6p2]

    cit3plot = []
    nio2plot = []
    nh4plot = []
    ni2pplot = []
    Hcitplot = []
    H2citplot = []
    H3citplot = []
    NiH2citplot = []
    NiHcitplot = []
    Nicitplot = []
    nh3plot = []
    nin4p2plot = []
    nin6p2plot = []
    for pHval in pH_x:
        cit3plot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[0])
        nio2plot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[1])
        nh4plot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[2])
        ni2pplot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[3])
        Hcitplot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[4])
        H2citplot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[5])
        H3citplot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[6])
        NiH2citplot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[7])
        NiHcitplot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[8])
        Nicitplot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[9])
        nh3plot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[10])
        nin4p2plot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[11])
        nin6p2plot.append(concs(citrate_total, ni_total, ammonia_total, pHval)[12])

    datasets= [cit3plot, nio2plot, nh4plot, ni2pplot, Hcitplot, H2citplot, H3citplot, NiH2citplot, NiHcitplot,
             Nicitplot,nh3plot, nin4p2plot, nin6p2plot]
    name = ['cit3', 'nio2', 'nh4', 'Ni2+', 'Hcit', 'H2cit', 'H3cit', 'NiH2cit', 'NiHcit', 'Nicit', 'nh3', 'Ni(Nh4)4 2+',
            'Ni(NH4)6 2+']
    fill = [None, None, None, None, None, None, None, None, None, None, None, None,None, None, None]
    color = ['rgb(90, 0, 100)', 'rgb(40, 130, 80)', 'rgb(245, 137, 22)', 'rgb(63, 63, 191)', 'rgb(191, 63, 63)', 'rgb(15, 15, 15)',
             'rgb(235, 64, 52)','rgb(137, 232, 186)','rgb(204, 75, 131)','rgb(14, 10, 209)','rgb(172, 51, 232)','rgb(2, 92, 8)','rgb(219, 140, 176)']

    data1 = []
    for i, dataset in enumerate(datasets):
        data1.append(go.Scatter(
            x=pH_x,
            y=dataset,
            mode='lines',
            hoverinfo='skip',
            fill=fill[i],
            name=name[i],
            showlegend=True,
            line=dict(
                shape='spline',
                width=2.5,
                color=color[i]
            )
        ))

    layout = go.Layout(
        xaxis={'title': 'pH', 'linecolor': 'grey', 'mirror':True},
        yaxis={'title': 'Concentration (kmolm<sup>-3</sup>)', 'linecolor': 'grey', 'mirror':True},
        # transition = {'duration': 1200},
        font=dict(family='Courier Sans', color='grey'),
        margin={'t': 50, 'l':10},
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgb(240,240,240)',
        autosize=False,
        width=600,
        height=500,
    )

    fig1 = go.Figure(data=data1, layout=layout)
    fig1.update_xaxes(gridcolor='LightPink', range=[0, 14],
                     nticks=20, mirror=True, ticks='outside', showline=True)

    fig1.update_yaxes(gridcolor='LightPink', ticks='outside',
                    range=[0, max(ni_total, citrate_total, ammonia_total,)*1.05])
    fig1.update_layout(
        title={
            'text': "Speciation plot",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig1
