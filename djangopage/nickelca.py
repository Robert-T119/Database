from djangopage.reactions3 import *
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
            value=1.2,
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
                value=1.3,
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
            value=1.4,
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

    datasets= [cit3plot, nio2plot, nh4plot, ni2pplot, Hcitplot, H2citplot, H3citplot, NiH2citplot, NiHcitplot,Nicitplot,nh3plot, nin4p2plot, nin6p2plot]
    name = ['cit3', 'nio2', 'nh4', 'Ni2+', 'Hcit', 'H2cit', 'H3cit', 'NiH2cit', 'NiHcit', 'Nicit', 'nh3', 'Ni(Nh4)4 2+','Ni(NH4)6 2+']
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

@app.callback(
    Output('Potential-pH curve', 'figure'),
    [Input('nickel_slider', 'value'),
     Input('citrate_slider', 'value'),
     Input('ammonia_dropdown', 'value')])

def speciation_graph(ni_total, citrate_total, ammonia_total):
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
    T_ = 298

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
     ###########################################################################################################################
    NiH2cit = NiHcit = Nicit = ni2pfree = nio2 = ni_total
    cit3 = Hcit = H2cit = H3cit = citrate_total
    nh3 = nh4 = ammonia_total
    if ni2pfree <= nh4 / 4:
        nin4p2 = ni2pfree
    else:
        nin4p2 = nh4 / 4
    if ni2pfree <= nh4 / 6:
        nin6p2 = ni2pfree
    else:
        nin6p2 = nh4 / 6
    T_ = 298
    pH_x = np.linspace(0, 14, 71)

    def nio2checker(nin4p2plot, nin6p2plot, ni2pplot, ammonia_total):
        if ammonia_total != 0.0:
            maxnin4p2 = max(nin4p2plot)
            i = nin4p2plot.index(maxnin4p2)
            maxnin6p2 = max(nin6p2plot)
            j = nin6p2plot.index(maxnin6p2)
            if maxnin6p2 > ni2pplot[j] and maxnin4p2 > ni2pplot[i]:
                status = 0
            elif maxnin6p2 > ni2pplot[j] and maxnin4p2 < ni2pplot[i]:
                status = 1
            elif maxnin6p2 < ni2pplot[j] and maxnin4p2 > ni2pplot[i]:
                status = 2
            elif maxnin6p2 < ni2pplot[j] and maxnin4p2 < ni2pplot[i]:
                status = 3
        elif ammonia_total == 0.0:
            status = 3
        return status
    status = nio2checker(nin4p2plot, nin6p2plot, ni2pplot, ammonia_total)
    # nio2checker(nin4p2plot, nin6p2plot, nip2pptplot)
    # nio2checker returns status. This can take 4 values, depending on what will be shown inside
    # the ni(oh)2 region
    # --------------------------------------------------------------------------------------------------------------
    def trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, nin4p2, nin6p2,nh4,nh3,nio2, T_):
        def interceptgenerator1(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_):
            bottom_eqs_1 = [[R1(ni2pfree, T_) for i in range(len(pH_x))], R6(pH_x, T_, NiH2cit, H3cit),R7(pH_x, T_, NiHcit, H3cit),
                            R8(pH_x, T_, NiHcit, H2cit), R9(pH_x, T_, Nicit, H2cit), R10(pH_x, T_, Nicit, Hcit),
                            [R11(T_, Nicit, cit3) for i in range(len(pH_x))], R12(pH_x, T_)]
            y_vars_1 = []
            # y vars has m, +c of every line in bottom eqs
            for i, eq in enumerate(bottom_eqs_1):
                y_vars_1.append(np.polyfit(pH_x, eq, 1))  # 1 means linear equations

            As_1 = []  # M
            Bs_1 = []  # C
            for i, var in enumerate(y_vars_1):
                if i >= 1:
                    As_1.append(np.array([[-y_vars_1[i - 1][0], 1], [-y_vars_1[i][0], 1]]))
                    Bs_1.append(np.array([y_vars_1[i - 1][1], y_vars_1[i][1]]))
            # inters has [x_intercept, y_intercept] of all intercepts
            inters_1 = []
            for i, ms in enumerate(As_1):
                for j, cs in enumerate(Bs_1):
                    if i == j:
                        # inters.append(np.linalg.solve(As[i],Bs[j]))
                        inters_1.append(np.linalg.inv(As_1[i]).dot(Bs_1[j]))
            return inters_1
        inters_1 = interceptgenerator1(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)

        def interceptgenerator2(pH_x, ni2pfree, nin4p2, nin6p2, nh4, nh3, T_):
            if status == 0:
                bottom_eqs_2 = [[R1(ni2pfree, T_) for i in range(len(pH_x))], R2(pH_x, T_), R3(pH_x, nh4, nin4p2, T_),
                                R4(pH_x, nh4, nin6p2, T_), [R5(nh3, nin6p2, T_) for i in range(len(pH_x))],R2(pH_x, T_)]
            elif status == 1:
                bottom_eqs_2 = [[R1(ni2pfree, T_) for i in range(len(pH_x))], R2(pH_x, T_),
                                R4(pH_x, nh4, nin6p2, T_), [R5(nh3, nin6p2, T_) for i in range(len(pH_x))],R2(pH_x, T_)]
            elif status == 2:
                bottom_eqs_2 = [[R1(ni2pfree, T_) for i in range(len(pH_x))], R2(pH_x, T_), R3(pH_x, nh4, nin4p2, T_)]

            elif status == 3 or ammonia_total == 0.0:
                bottom_eqs_2 = [[R1(ni2pfree, T_) for i in range(len(pH_x))], R2(pH_x, T_)]

            y_vars_2 = []
            # y vars has m, +c of every line in bottom eqs
            for i, eq in enumerate(bottom_eqs_2):
                y_vars_2.append(np.polyfit(pH_x, eq, 1))

            As_2 = []
            Bs_2 = []
            for i, var in enumerate(y_vars_2):
                if i >= 1:
                    As_2.append(np.array([[-y_vars_2[i - 1][0], 1], [-y_vars_2[i][0], 1]]))
                    Bs_2.append(np.array([y_vars_2[i - 1][1], y_vars_2[i][1]]))
            # inters has [x_intercept, y_intercept] of all intercepts
            inters_2 = []
            for i, ms in enumerate(As_2):
                for j, cs in enumerate(Bs_2):
                    if i == j:
                        inters_2.append(np.linalg.inv(As_2[i]).dot(Bs_2[j]))
            return inters_2
        inters_2 = interceptgenerator2(pH_x, ni2pfree, nin4p2, nin6p2, nh4, nh3, T_)

        xinters_1 = []
        for item in inters_1:
            xinters_1.append(item[0])

        x_data_1 = []
        for i, item in enumerate(xinters_1):
            if i == 0:
                x_data_1.append(np.linspace(0, item, 5))
            elif i >= 1:
                x_data_1.append(np.linspace(xinters_1[i - 1], item, 5))
        finalindex_1 = len(xinters_1) - 1
        x_data_1.append(np.linspace(xinters_1[finalindex_1], 14, 5))

        xinters_2 = []
        for item in inters_2:
            xinters_2.append(item[0])
        x_data_2 = []
        if status != 3:
            for i, item in enumerate(xinters_2):
                if i == 0:
                    x_data_2.append(np.linspace(0, item, 5))
                elif i >= 1:
                    x_data_2.append(np.linspace(xinters_2[i - 1], item, 5))
            finalindex_2 = len(xinters_2) - 1
            x_data_2.append(np.linspace(xinters_2[finalindex_2], 14, 5))
        elif status == 3 or ammonia_total == 0.0:
            x_data_2.append(list(np.linspace(0, xinters_2[0], 5)))
            x_data_2.append(list(np.linspace(xinters_2[0], 14, 5)))

        y_data_bottom_1 = [[R1(ni2pfree, T_) for i in range(len(x_data_1[0]))], R6(x_data_1[1], T_, NiH2cit, H3cit),
                           R7(x_data_1[2], T_, NiHcit, H3cit),
                           R8(x_data_1[3], T_, NiHcit, H2cit), R9(x_data_1[4], T_, Nicit, H2cit),
                           R10(x_data_1[5], T_, Nicit, Hcit),
                           [R11(T_, Nicit, cit3) for i in range(len(x_data_1[6]))], R12(x_data_1[7], T_)]

        new_x_bottom_1 = []
        new_y_bottom_1 = []

        for xvalues_1 in x_data_1:
            for xvalue_1 in xvalues_1:
                new_x_bottom_1.append(xvalue_1)
        for yvalues_1 in y_data_bottom_1:
            for yvalue_1 in yvalues_1:
                new_y_bottom_1.append(yvalue_1)

        if status == 0:
            y_data_bottom_2 = [[R1(ni2pfree, T_) for i in range(len(x_data_2[0]))], R2(x_data_2[1], T_),
                     R3(x_data_2[2], nh4, nin4p2, T_),R4(x_data_2[3], nh4, nin6p2, T_), [R5(nh3, nin6p2, T_) for i in range(len(x_data_2[4]))],R2(x_data_2[5], T_)]
        elif status == 1:
            y_data_bottom_2 = [[R1(ni2pfree, T_) for i in range(len(x_data_2[0]))], R2(x_data_2[1], T_),
                               R4(x_data_2[2], nh4, nin6p2, T_), [R5(nh3, nin6p2, T_) for i in range(len(x_data_2[3]))],R2(x_data_2[4], T_)]
        elif status == 2:
            y_data_bottom_2 = [[R1(ni2pfree, T_) for i in range(len(x_data_2[0]))], R2(x_data_2[1], T_), R3(x_data_2[2], nh4, nin4p2, T_),R2(x_data_2[3], T_)]
        elif status == 3 or ammonia_total == 0.0:
            y_data_bottom_2 = [[R1(ni2pfree, T_) for i in range(len(x_data_2[0]))], R2(x_data_2[1], T_)]

        new_x_bottom_2 = []
        new_y_bottom_2 = []
        for xvalues in x_data_2:
            for xvalue in xvalues:
                new_x_bottom_2.append(xvalue)
        for yvalues in y_data_bottom_2:
            for yvalue in yvalues:
                new_y_bottom_2.append(yvalue)

        y_data_top_1 = [T1(ni2pfree, x_data_1[0], T_), T6(x_data_1[1], NiH2cit, H3cit, T_),
                        T7(x_data_1[2], NiHcit, H3cit, T_),
                        T8(x_data_1[3], NiHcit, H2cit, T_),
                        T9(x_data_1[4], Nicit, H2cit, T_), T10(x_data_1[5], Nicit, Hcit, T_),
                        T11(x_data_1[6], Nicit, cit3, T_),
                        T12(x_data_1[7], T_)]
        new_x_top_1 = []
        new_y_top_1 = []
        for xvalues in x_data_1:
            for xvalue in xvalues:
                new_x_top_1.append(xvalue)
        for yvalues in y_data_top_1:
            for yvalue in yvalues:
                new_y_top_1.append(yvalue)

        if status == 0:
            y_data_top_2 = [T1(ni2pfree, x_data_2[0], T_), T2(x_data_2[1], T_), T3(nh4, nin4p2, x_data_2[2], T_),
                            T4(x_data_2[3], nh4, nin6p2, T_), T5(nh4, x_data_2[4], nin6p2, T_), T2(x_data_2[5], T_)]
        elif status == 1:
            y_data_top_2 = [T1(ni2pfree, x_data_2[0], T_), T2(x_data_2[1], T_),
                            T4(x_data_2[2], nh4, nin6p2, T_), T5(nh4, x_data_2[3], nin6p2, T_), T2(x_data_2[4], T_)]
        elif status == 2:
            y_data_top_2 = [T1(ni2pfree, x_data_2[0], T_), T2(x_data_2[1], T_), T3(nh4, nin4p2, x_data_2[2], T_),
                          T2(x_data_2[3], T_)]
        elif status == 3 or ammonia_total == 0.0:
            y_data_top_2 = [T1(ni2pfree, x_data_2[0], T_), T2(x_data_2[1], T_)]

        new_x_top_2 = []
        new_y_top_2 = []
        for xvalues in x_data_2:
            for xvalue in xvalues:
                new_x_top_2.append(xvalue)
        for yvalues in y_data_top_2:
            for yvalue in yvalues:
                new_y_top_2.append(yvalue)

        y_interps_1 = [T1(ni2pfree, inters_1[0][0], T_), T6(inters_1[1][0], NiH2cit, H3cit, T_),
                       T7(inters_1[2][0], NiHcit, H3cit, T_), T8(inters_1[3][0], NiHcit, H2cit, T_),
                       T9(inters_1[4][0], Nicit, H2cit, T_), T10(inters_1[5][0], Nicit, Hcit, T_),
                       T11(inters_1[6][0], Nicit, cit3, T_)]
        vys_1 = []
        for i, val in enumerate(inters_1):
            vys_1.append(list(np.linspace(inters_1[i][1], y_interps_1[i], 5)))
        x_data_verticals_1 = []
        for i in range(len(inters_1)):
            x_data_verticals_1.append([inters_1[i][0] for j in range(len(vys_1[i]))])
        new_x_vert_1 = []
        new_y_vert_1 = []
        for xvalues in x_data_verticals_1:
            for xvalue in xvalues:
                new_x_vert_1.append(xvalue)
        for yvalues in vys_1:
            for yvalue in yvalues:
                new_y_vert_1.append(yvalue)

        if status == 0:
            y_interps_2 = [T1(ni2pfree, inters_2[0][0], T_), T2(inters_2[1][0], T_),
                           T3(nh4, nin4p2, inters_2[2][0], T_),T4(inters_2[3][0], nh4, nin6p2, T_), T5(nh4, inters_2[4][0], nin6p2, T_)]
        elif status == 1:
            y_interps_2 = [T1(ni2pfree, inters_2[0][0], T_), T2(inters_2[1][0], T_),T4(inters_2[2][0], nh4, nin6p2, T_), T5(nh4, inters_2[3][0], nin6p2, T_)]
        elif status == 2:
            y_interps_2 = [T1(ni2pfree, inters_2[0][0], T_), T2(inters_2[1][0], T_), T3(nh4, nin4p2, inters_2[2][0], T_)]
        elif status == 3 or ammonia_total == 0.0:
            y_interps_2 = [T1(ni2pfree, inters_2[0][0], T_)]
        vys_2 = []
        for i, val in enumerate(inters_2):
            vys_2.append(list(np.linspace(inters_2[i][1], y_interps_2[i], 5)))
        x_data_verticals_2 = []
        for i in range(len(inters_2)):
            x_data_verticals_2.append([inters_2[i][0] for j in range(len(vys_2[i]))])
        new_x_vert_2 = []
        new_y_vert_2 = []
        for xvalues in x_data_verticals_2:
            for xvalue in xvalues:
                new_x_vert_2.append(xvalue)
        for yvalues in vys_2:
            for yvalue in yvalues:
                new_y_vert_2.append(yvalue)

        if status == 0:
            nio3regionx = list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list(
                [0 for i in range(0, 5)]) + list(x_data_1[0]) + list(x_data_1[1]) + list(x_data_1[2]) + list(
                np.linspace(inters_1[5][0], inters_2[2][0], 5)) + list(x_data_2[3]) + list(x_data_2[4]) + list(
                x_data_2[5])
            nio3regiony = list(np.linspace(T12(14, T_), 2.6, 5)) + list([2.6 for i in range(0, 5)]) + list(
                np.linspace(T1(ni2pfree, 0, 298), 2.6, 5)) + list(y_data_top_1[0]) + list(y_data_top_1[1]) + list(
                y_data_top_1[2]) + list(
                np.linspace(T11(inters_1[5][0], Nicit, cit3, T_), T4(inters_2[2][0], nh4, nin6p2, T_), 5)) + list(
                y_data_top_2[3]) + list(y_data_top_2[4]) + list(y_data_top_2[5])

            niregionx = list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list(
                [0 for i in range(0, 5)]) + list(x_data_1[0]) + list(x_data_1[1]) + list(x_data_1[2]) + list(
                np.linspace(inters_1[5][0], inters_2[2][0], 5)) + list(x_data_2[3]) + list(x_data_2[4]) + list(
                x_data_2[5])
            niregiony = list(np.linspace(R12(14, T_), -1.8, 5)) + list([-1.8 for i in range(0, 5)]) + list(
                np.linspace(-1.8, R1(ni2pfree, 298), 5)) + list(y_data_bottom_1[0]) + list(y_data_bottom_1[1]) + list(
                y_data_bottom_1[2]) + list(np.linspace(inters_1[5][1], inters_2[2][1], 5)) + list(
                y_data_bottom_2[3]) + list(y_data_bottom_2[4]) + list(y_data_bottom_2[5])

            nip2regionx = list(x_data_1[0]) + list(x_data_verticals_1[0]) + list(reversed(x_data_1[0])) + list(
                [0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom_1[0]) + list(vys_1[0]) + list(reversed(y_data_top_1[0])) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, 0, 298), 5))

            NiH2citregionx = list(x_data_1[1]) + list(x_data_1[2]) + list(x_data_verticals_1[0]) + list(
                x_data_1[1]) + list(
                x_data_1[2]) + list(x_data_verticals_1[1])
            NiH2citregiony = list(y_data_bottom_1[1]) + list(y_data_bottom_1[2]) + list(vys_1[0]) + list(
                y_data_top_1[1]) + list(
                y_data_top_1[2]) + list(reversed(vys_1[1]))

            NiHcitregionx = list(x_data_1[3]) + list(x_data_verticals_1[1]) + list(x_data_1[3]) + list(
                x_data_verticals_1[3])
            NiHcitregiony = list(y_data_bottom_1[3]) + list(vys_1[1]) + list(y_data_top_1[3]) + list(reversed(vys_1[3]))

            Nicitregionx = list(x_data_verticals_1[3]) + list(x_data_1[4]) + list(x_data_1[5]) + list(
                np.linspace(inters_1[5][0], inters_2[2][0], 5)) + list(reversed(x_data_2[1])) + list(
                x_data_verticals_2[0]) + list(x_data_2[1]) + list(
                reversed(np.linspace(inters_1[5][0], inters_2[2][0], 5))) + list(reversed(x_data_1[5])) + list(
                reversed(x_data_1[4]))
            Nicitregiony = list(vys_1[3]) + list(y_data_bottom_1[4]) + list(y_data_bottom_1[5]) + list(
                np.linspace(inters_1[5][1], inters_2[2][1], 5)) + list(reversed(y_data_bottom_2[1])) + list(
                vys_2[0]) + list(y_data_top_2[1]) + list(reversed(
                np.linspace(T11(inters_1[5][0], Nicit, cit3, T_), T4(inters_2[2][0], nh4, nin6p2, T_), 5))) + list(
                reversed(y_data_top_1[5])) + list(reversed(y_data_top_1[4]))

            nio2regionx1 = list(x_data_2[1]) + list(x_data_verticals_2[0]) + list(x_data_2[1]) + list(
                x_data_verticals_2[1])
            nio2regiony1 = list(y_data_bottom_2[1]) + list(vys_2[0]) + list(y_data_top_2[1]) + list(reversed(vys_2[1]))

            nin4p2regionx = list(x_data_verticals_2[1]) + list(x_data_2[2]) + list(x_data_verticals_2[2]) + list(
                reversed(x_data_2[2]))
            nin4p2regiony = list(vys_2[2]) + list(y_data_bottom_2[2]) + list(reversed(vys_2[2])) + list(y_data_top_2[2])

            nin6p2regionx = list(x_data_verticals_2[2]) + list(x_data_2[3]) + list(x_data_2[4]) + list(
                x_data_verticals_2[4]) + list(reversed(x_data_2[4])) + list(reversed(x_data_2[3]))
            nin6p2regiony = list(vys_2[2]) + list(y_data_top_2[3]) + list(y_data_top_2[4]) + list(
                reversed(vys_2[4])) + list(reversed(y_data_bottom_2[4])) + list(reversed(y_data_bottom_2[3]))

            nio2regionx2 = list(reversed(x_data_2[5])) + list(x_data_verticals_2[4]) + list(x_data_2[5]) + list(
                14 for i in range(0, 5))
            nio2regiony2 = list(reversed(y_data_bottom_2[5])) + list(vys_2[4]) + list(y_data_top_2[5]) + list(
                np.linspace(R12(14, 298), T12(14, 298), 5))

            xs = [niregionx,nio3regionx,nip2regionx, nio2regionx1, NiH2citregionx, NiHcitregionx, Nicitregionx, nin4p2regionx, nin6p2regionx,
                  nio2regionx2]
            ys = [niregiony,nio3regiony,nip2regiony, nio2regiony1, NiH2citregiony, NiHcitregiony, Nicitregiony, nin4p2regiony, nin6p2regiony,
                  nio2regiony2]

        elif status == 1:
            nio3regionx = list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list(
                [0 for i in range(0, 5)]) + list(x_data_1[0]) + list(x_data_1[1]) + list(x_data_1[2]) + list(
                np.linspace(inters_1[5][0], inters_2[2][0], 5)) + list(x_data_2[3]) + list(x_data_2[4]) + list(
                x_data_2[5])
            nio3regiony = list(np.linspace(T12(14, T_), 2.6, 5)) + list([2.6 for i in range(0, 5)]) + list(
                np.linspace(T1(ni2pfree, 0, 298), 2.6, 5)) + list(y_data_top_1[0]) + list(y_data_top_1[1]) + list(
                y_data_top_1[2]) + list(
                np.linspace(T11(inters_1[5][0], Nicit, cit3, T_), T4(inters_2[2][0], nh4, nin6p2, T_), 5)) + list(
                y_data_top_2[3]) + list(y_data_top_2[4]) + list(y_data_top_2[5])

            niregionx = list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list(
                [0 for i in range(0, 5)]) + list(x_data_1[0]) + list(x_data_1[1]) + list(x_data_1[2]) + list(
                np.linspace(inters_1[5][0], inters_2[2][0], 5)) + list(x_data_2[3]) + list(x_data_2[4]) + list(
                x_data_2[5])
            niregiony = list(np.linspace(R12(14, T_), -1.8, 5)) + list([-1.8 for i in range(0, 5)]) + list(
                np.linspace(-1.8, R1(ni2pfree, 298), 5)) + list(y_data_bottom_1[0]) + list(y_data_bottom_1[1]) + list(
                y_data_bottom_1[2]) + list(np.linspace(inters_1[5][1], inters_2[2][1], 5)) + list(
                y_data_bottom_2[3]) + list(y_data_bottom_2[4]) + list(y_data_bottom_2[5])

            nip2regionx = list(x_data_1[0]) + list(x_data_verticals_1[0]) + list(reversed(x_data_1[0])) + list(
                [0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom_1[0]) + list(vys_1[0]) + list(reversed(y_data_top_1[0])) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, 0, 298), 5))

            NiH2citregionx = list(x_data_1[1]) + list(x_data_1[2]) + list(x_data_verticals_1[0]) + list(
                x_data_1[1]) + list(
                x_data_1[2]) + list(x_data_verticals_1[1])
            NiH2citregiony = list(y_data_bottom_1[1]) + list(y_data_bottom_1[2]) + list(vys_1[0]) + list(
                y_data_top_1[1]) + list(
                y_data_top_1[2]) + list(reversed(vys_1[1]))

            NiHcitregionx = list(x_data_1[3]) + list(x_data_verticals_1[1]) + list(x_data_1[3]) + list(
                x_data_verticals_1[3])
            NiHcitregiony = list(y_data_bottom_1[3]) + list(vys_1[1]) + list(y_data_top_1[3]) + list(reversed(vys_1[3]))

            Nicitregionx = list(x_data_verticals_1[3]) + list(x_data_1[4]) + list(x_data_1[5]) + list(
                np.linspace(inters_1[5][0], inters_2[2][0], 5)) + list(reversed(x_data_2[1])) + list(
                x_data_verticals_2[0]) + list(x_data_2[1]) + list(
                reversed(np.linspace(inters_1[5][0], inters_2[2][0], 5))) + list(reversed(x_data_1[5])) + list(
                reversed(x_data_1[4]))
            Nicitregiony = list(vys_1[3]) + list(y_data_bottom_1[4]) + list(y_data_bottom_1[5]) + list(
                np.linspace(inters_1[5][1], inters_2[2][1], 5)) + list(reversed(y_data_bottom_2[1])) + list(
                vys_2[0]) + list(y_data_top_2[1]) + list(reversed(
                np.linspace(T11(inters_1[5][0], Nicit, cit3, T_), T4(inters_2[2][0], nh4, nin6p2, T_), 5))) + list(
                reversed(y_data_top_1[5])) + list(reversed(y_data_top_1[4]))

            nio2regionx1 = list(x_data_2[1]) + list(x_data_verticals_2[0]) + list(x_data_2[1]) + list(
                x_data_verticals_2[1])
            nio2regiony1 = list(y_data_bottom_2[1]) + list(vys_2[0]) + list(y_data_top_2[1]) + list(reversed(vys_2[1]))

            nin6p2regionx = list(x_data_verticals_2[1]) + list(x_data_2[2]) + list(x_data_2[3]) + list(
                x_data_verticals_2[3]) + list(reversed(x_data_2[3])) + list(reversed(x_data_2[2]))
            nin6p2regiony = list(vys_2[1]) + list(y_data_top_2[2]) + list(y_data_top_2[3]) + list(
                reversed(vys_2[3])) + list(reversed(y_data_bottom_2[3])) + list(reversed(y_data_bottom_2[2]))

            nio2regionx2 = list(reversed(x_data_2[4])) + list(x_data_verticals_2[3]) + list(x_data_2[4]) + list(
                14 for i in range(0, 5))
            nio2regiony2 = list(reversed(y_data_bottom_2[4])) + list(vys_2[3]) + list(y_data_top_2[4]) + list(
                np.linspace(R12(14, 298), T12(14, 298), 5))

            xs = [niregionx,nio3regionx,nip2regionx,nio2regionx1, NiH2citregionx, NiHcitregionx, Nicitregionx, nin6p2regionx, nio2regionx2]
            ys = [niregiony,nio3regiony,nip2regiony,nio2regiony1, NiH2citregiony, NiHcitregiony, Nicitregiony, nin6p2regiony, nio2regiony2]

        elif status == 2:
            nio3regionx = list(new_x_bottom_1) + list([14 for i in range(0, 5)]) + list(
                reversed(np.linspace(0, 14, 5))) + list(
                [0 for i in range(0, 5)])
            nio3regiony = list(new_y_top_1) + list(np.linspace(T12(14, T_), 2.6, 5)) + list(
                [2.6 for i in range(0, 5)]) + list(
                np.linspace(T1(ni2pfree, 0, 298), 2.6, 5))

            niregionx = list(new_x_bottom_1) + list([14 for i in range(0, 5)]) + list(
                reversed(np.linspace(0, 14, 5))) + list(
                [0 for i in range(0, 5)])
            niregiony = list(new_y_bottom_1) + list(np.linspace(R12(14, T_), -1.8, 5)) + list(
                [-1.8 for i in range(0, 5)]) + list(
                np.linspace(-1.8, R1(ni2pfree, 298), 5))

            nip2regionx = list(x_data_1[0]) + list(x_data_verticals_1[0]) + list(reversed(x_data_1[0])) + list(
                [0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom_1[0]) + list(vys_1[0]) + list(reversed(y_data_top_1[0])) + list(
                np.linspace(R1(nip2regionx, T_), T1(ni2pfree, 0, 298), 5))

            NiH2citregionx = list(x_data_1[1]) + list(x_data_1[2]) + list(x_data_verticals_1[0]) + list(x_data_1[1]) + list(
                x_data_1[2]) + list(x_data_verticals_1[1])
            NiH2citregiony = list(y_data_bottom_1[1]) + list(y_data_bottom_1[2]) + list(vys_1[0]) + list(
                y_data_top_1[1]) + list(y_data_top_1[2]) + list(reversed(vys_1[1]))

            NiHcitregionx = list(x_data_1[3]) + list(x_data_verticals_1[1]) + list(x_data_1[3]) + list(x_data_verticals_1[3])
            NiHcitregiony = list(y_data_bottom_1[3]) + list(vys_1[1]) + list(y_data_top_1[3]) + list(reversed(vys_1[3]))

            Nicitregionx = list(x_data_1[4]) + list(x_data_1[5]) + list(x_data_1[6]) + list(x_data_verticals_1[3]) + list(
                x_data_1[4]) + list(x_data_1[5]) + list(x_data_1[6]) + list(x_data_verticals_1[6])
            Nicitregiony = list(y_data_bottom_1[4]) + list(y_data_bottom_1[5]) + list(y_data_bottom_1[6]) + list(
                vys_1[3]) + list(y_data_top_1[4]) + list(y_data_top_1[5]) + list(y_data_top_1[6]) + list(reversed(vys_1[6]))

            nio2regionx = list(reversed(x_data_1[7])) + list(x_data_verticals_1[6]) + list(x_data_1[7]) + list(
                [14 for i in range(0, 5)])
            nio2regiony = list(reversed(y_data_bottom_1[7])) + list(vys_1[6]) + list(y_data_top_1[7]) + list(
                np.linspace(R12(14, T_), T12(14, T_), 5))
            xs = [niregionx,nio3regiony,nip2regionx, NiH2citregionx, NiHcitregionx, Nicitregionx, nio2regionx]
            ys = [niregionx,nio3regiony,nip2regiony, NiH2citregiony, NiHcitregiony, Nicitregiony, nio2regiony]

        elif status == 3:
            nio3regionx = list(new_x_bottom_1) + list([14 for i in range(0, 5)]) + list(
                reversed(np.linspace(0, 14, 5))) + list(
                [0 for i in range(0, 5)])
            nio3regiony = list(new_y_top_1) + list(np.linspace(T12(14, T_), 2.6, 5)) + list(
                [2.6 for i in range(0, 5)]) + list(
                np.linspace(T1(ni2pfree, 0, 298), 2.6, 5))

            niregionx = list(new_x_bottom_1) + list([14 for i in range(0, 5)]) + list(
                reversed(np.linspace(0, 14, 5))) + list(
                [0 for i in range(0, 5)])
            niregiony = list(new_y_bottom_1) + list(np.linspace(R12(14, T_), -1.8, 5)) + list(
                [-1.8 for i in range(0, 5)]) + list(
                np.linspace(-1.8, R1(ni2pfree, 298), 5))

            nip2regionx = list(x_data_1[0]) + list(x_data_verticals_1[0]) + list(reversed(x_data_1[0])) + list(
                [0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom_1[0]) + list(vys_1[0]) + list(reversed(y_data_top_1[0])) + list(
                np.linspace(R1(nip2regionx, T_), T1(ni2pfree, 0, 298), 5))

            NiH2citregionx = list(x_data_1[1]) + list(x_data_1[2]) + list(x_data_verticals_1[0]) + list(
                x_data_1[1]) + list(
                x_data_1[2]) + list(x_data_verticals_1[1])
            NiH2citregiony = list(y_data_bottom_1[1]) + list(y_data_bottom_1[2]) + list(vys_1[0]) + list(
                y_data_top_1[1]) + list(y_data_top_1[2]) + list(reversed(vys_1[1]))

            NiHcitregionx = list(x_data_1[3]) + list(x_data_verticals_1[1]) + list(x_data_1[3]) + list(
                x_data_verticals_1[3])
            NiHcitregiony = list(y_data_bottom_1[3]) + list(vys_1[1]) + list(y_data_top_1[3]) + list(reversed(vys_1[3]))

            Nicitregionx = list(x_data_1[4]) + list(x_data_1[5]) + list(x_data_1[6]) + list(
                x_data_verticals_1[3]) + list(
                x_data_1[4]) + list(x_data_1[5]) + list(x_data_1[6]) + list(x_data_verticals_1[6])
            Nicitregiony = list(y_data_bottom_1[4]) + list(y_data_bottom_1[5]) + list(y_data_bottom_1[6]) + list(
                vys_1[3]) + list(y_data_top_1[4]) + list(y_data_top_1[5]) + list(y_data_top_1[6]) + list(
                reversed(vys_1[6]))

            nio2regionx = list(reversed(x_data_1[7])) + list(x_data_verticals_1[6]) + list(x_data_1[7]) + list(
                [14 for i in range(0, 5)])
            nio2regiony = list(reversed(y_data_bottom_1[7])) + list(vys_1[6]) + list(y_data_top_1[7]) + list(
                np.linspace(R12(14, T_), T12(14, T_), 5))
            xs = [niregionx,nio3regiony,nip2regionx, NiH2citregionx, NiHcitregionx, Nicitregionx, nio2regionx]
            ys = [niregionx,nio3regiony,nip2regiony, NiH2citregiony, NiHcitregiony, Nicitregiony, nio2regiony]
        return [xs, ys]
    # end function, return data to add to traces, should already be in correct form to
    # cooperate with dash notation
    #------------------------------------------------------------------------------------------------
    xs =trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, nin4p2, nin6p2, nh4, nh3, nio2,T_)[0]
    ys =trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, nin4p2, nin6p2, nh4, nh3, nio2,T_)[1]

    if status == 0:
        name = ['Ni', 'Ni(OH)<sub>3</sub>','Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>', 'NiH<sub>2</sub>cit<sup>+</sup>', 'NiHcit</sub>',
                'Nicit<sup>-</sup>', '[Ni(NH<sub>3</sub>)<sub>4</sub>]<sup>2+</sup>',
                '[Ni(NH<sub>3</sub>)<sub>6</sub>]<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(127, 63, 191, 0.5)', 'rgba(30, 205, 40, 0.5)','rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)', 'rgba(114, 102, 234, 0.63)',
                 'rgba(114, 204, 234, 0.63)', 'rgba(245, 40, 145, 0.8)', 'rgba(87, 131, 200, 0.8)',
                 'rgba(8, 217, 167, 0.66)', 'rgba(243, 89, 83, 0.32)', 'rgba(222, 76, 63, 0.63)',
                 'rgba(55, 225, 234, 0.63)']

    elif status == 1:
        name = ['Ni', 'Ni(OH)<sub>3</sub>','Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>', 'NiH<sub>2</sub>cit<sup>+</sup>', 'NiHcit</sub>',
                'Nicit<sup>-</sup>','[Ni(NH<sub>3</sub>)<sub>6</sub>]<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(127, 63, 191, 0.5)', 'rgba(30, 205, 40, 0.5)','rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)', 'rgba(114, 102, 234, 0.63)',
                 'rgba(114, 204, 234, 0.63)', 'rgba(245, 40, 145, 0.8)', 'rgba(87, 131, 200, 0.8)',
                 'rgba(8, 217, 167, 0.66)', 'rgba(222, 76, 63, 0.63)','rgba(55, 225, 234, 0.63)']

    elif status == 2: # this never really happens
        name = ['Ni', 'Ni(OH)<sub>3</sub>','Ni<sup>2+</sup>', 'NiH<sub>2</sub>cit<sup>+</sup>', 'NiHcit</sub>', 'Nicit<sup>-</sup>','Ni(OH)<sub>2</sub>']

    elif status == 3:
        name = ['Ni', 'Ni(OH)<sub>3</sub>','Ni<sup>2+</sup>', 'NiH<sub>2</sub>cit<sup>+</sup>', 'NiHcit</sub>', 'Nicit<sup>-</sup>','Ni(OH)<sub>2</sub>']
        color = ['rgba(127, 63, 191, 0.5)', 'rgba(30, 205, 40, 0.5)','rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)', 'rgba(252, 177, 101, 0.8)',
                   'rgba(7, 117, 189, 0.66)', 'rgba(63, 63, 191, 0.5)']

    data = []
    if status != 3:
        for i, xvals in enumerate(xs):
            data.append(go.Scatter(
                x=xvals,
                y=ys[i],
                mode='none',
                fill='toself',
                hoverinfo='skip',
                fillcolor=color[i],
                showlegend=True,
                name=name[i]
            ))

    elif status == 3:
        for i, xvals in enumerate(xs):
            data.append(go.Scatter(
                x=xvals,
                y=ys[i],
                mode='none',
                fill='toself',
                hoverinfo='skip',
                fillcolor=color[i],
                showlegend=True,
                name=name[i]
            ))
    # add water splitting
    ywater = [W1(pH_x, T_), W2(pH_x, T_)]
    for ys in ywater:
        data.append(go.Scatter(
            x = pH_x,
            y = ys,
            mode='lines',
            hoverinfo='skip',
            showlegend=False,
            line=dict(
                shape='spline',
                color='blue',
                width=2,
                dash='dash'
                )
            ))

    layout = go.Layout(
        xaxis={'title': 'pH', 'linecolor':'grey', 'mirror':True},
        yaxis={'title': "Electrode potential, <i>E</i> vs SHE/ V ", 'linecolor':'grey', 'mirror':True},
        transition = {'duration': 1200},
        font=dict(family='Courier Sans', color='grey'),
        margin={'t': 50, 'r': 20},
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(50,0,0,0)',
        autosize=False,
        width=600,
        height=500,
            )

    fig = go.Figure(data=data, layout=layout)
    if ammonia_total != 0:
        fig.add_trace(go.Scatter(
            x=list([9.25 for i in range(len(np.linspace(-1, 2.4, 5)))]),
            y=list(np.linspace(-1, 2.4, 5)),
            mode='lines',
            showlegend=False,
            hoverinfo='skip',
            line=dict(
                    shape='spline',
                    color='purple',
                    width=1.5,
                    dash='dash'
                    )
        ))
        fig.add_trace(go.Scatter(
            x=list([2.996 for i in range(len(np.linspace(-1, 2.4, 5)))]),
            y=list(np.linspace(-1, 2.4, 5)),
            mode='lines',
            showlegend=False,
            hoverinfo='skip',
            line=dict(
                shape='spline',
                color='green',
                width=1.5,
                dash='dash'
            )
        ))
        fig.add_trace(go.Scatter(
            x=list([4.7478 for i in range(len(np.linspace(-1, 2.4, 5)))]),
            y=list(np.linspace(-1, 2.4, 5)),
            mode='lines',
            showlegend=False,
            hoverinfo='skip',
            line=dict(
                shape='spline',
                color='green',
                width=1.5,
                dash='dash'
            )
        ))
        fig.add_trace(go.Scatter(
            x=list([6.3946 for i in range(len(np.linspace(-1, 2.4, 5)))]),
            y=list(np.linspace(-1, 2.4, 5)),
            mode='lines',
            showlegend=False,
            hoverinfo='skip',
            line=dict(
                shape='spline',
                color='green',
                width=1.5,
                dash='dash'
            )
        ))
        fig.add_annotation(
            x=8.2,
            y=2,
            text="NH<sub>4</sub><sup>+</sup>",
            showarrow=False,
            font=dict(
                    family="Courier New bold",
                    size=20,
                    color="purple"
                    ),)
        fig.add_annotation(
            x=10.3,
            y=1.97,
            text="NH<sub>3</sub>",
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=20,
                color="purple"
            )
        )
        fig.add_annotation(
            x=2.1,
            y=2.2,
            text='H<sub>3</sub>cit',
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=16,
                color="purple"
            ), )
        fig.add_annotation(
            x=3.8,
            y=2.2,
            text='H<sub>2</sub>cit<sup>2-</sup>',
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=16,
                color="purple"
            )
        )
        fig.add_annotation(
            x=5.6,
            y=2.2,
            text='Hcit<sup>2-</sup>',
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=16,
                color="purple"
            )
        )
        fig.add_annotation(
            x=7.2,
            y=2.2,
            text='Cit<sup>3-</sup>',
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=16,
                color="purple"
            )
        )
    fig.add_annotation(
        x=0.4,
        y=1.37,
        showarrow=False,
        text='O<sub>2</sub>',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue',
        )
    )
    fig.add_annotation(
        x=0.55,
        y=1.02,
        showarrow=False,
        text='H<sub>2</sub>O',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue'
        )
    )
    fig.add_annotation(
        x=0.37,
        y=0.1,
        showarrow=False,
        text='H<sup>+</sup>',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue'
        )
    )
    fig.add_annotation(
        x=0.4,
        y=-0.2,
        showarrow=False,
        text='H<sub>2</sub>',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue'
        )
    )
    fig.update_xaxes(gridcolor='LightPink',
                     nticks=20,  ticks='outside', showline=True,
                    zeroline=True, zerolinewidth=2, zerolinecolor='LightPink')
    fig.update_yaxes(nticks=20, gridcolor='LightPink',  ticks='outside', showline=True,
                    zeroline=True, zerolinewidth=2, zerolinecolor='LightPink')
    fig.update_layout(yaxis=dict(range=[-1,2.4]), xaxis=dict(range=[0,14]),
                        title={
                            'text': "Potential-pH diagram",
                            'y':0.95,
                            'x':0.5,
                            'xanchor': 'center',
                            'yanchor': 'top'
                            })
    return fig