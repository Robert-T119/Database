from djangopage.reactions2 import *
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import fsolve
from django_plotly_dash import DjangoDash

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app = DjangoDash('nickelcitrate', add_bootstrap_links=True)
# , add_bootstrap_links=True
# add_bootstrap_links=True


app.layout = html.Div([
    html.Div([
        html.H2(u'Nickel-Citrate-H\u2082O system')],
        style={
            'text-align':'center',
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '-10px 0 -10px 0',
            'margin-bottom': '2px',
        }
        ),

    html.Div([
        html.H6(u"Total nickel concentration (kmolm\u207B\u00B3):"),
        dcc.Slider(
            id='nickel_slider',
            min=0.1,
            max=3.0,
            value=1.0,
            step=0.1,
            marks={n_activity: str(n_activity) for n_activity in [0.1,0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3]},
                ),
            ],
        style={
            'padding': '10px 15px 10px',
            'border': 'thin lightgrey solid',
            'margin-bottom': '3px',
            'backgroundColor': 'rgb(250, 250, 250)',
            }
        ),

    html.Div([
        html.H6(u'Total citrate concentration (kmolm\u207B\u00B3):'),
        dcc.Slider(
            id='citrate_dropdown',
            min=0.1,
            max=3.0,
            value=1.0,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1,0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3]},
                 ),
    ],
        style={
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250, 250, 250)',
            'padding': '10px 15px 10px',
            'margin-bottom': '3px',
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
    Output('Potential-pH curve', 'figure'),
    [Input('nickel_slider', 'value'),
     Input('citrate_dropdown', 'value')])

def speciation_graph(citrate_total, ni_total):
    k1 = 9.197 * (10 ** 11)
    k2 = 1.01 * (10 ** -3)
    k3 = 1.787 * (10 ** -5)
    k4 = 4.031 * (10 ** -7)
    k5 = 2.512 * (10 ** 5)
    k6 = 1.995 * (10 ** 3)
    k7 = 5.623 * (10 ** 1)
    logk = 11.96
    pH_x = np.linspace(0, 14, 71)
    #----------------------------------------------------------------------------------------------
    # begin first function, output all species concentrations. One concentration for each pH value.

    def concs(citrate_total, ni_total, pH_x):
        h = 10 ** (-pH_x)

        if citrate_total != 0:
            def f(z):
                cit3 = z[0]
                nio2 = z[1]
                F = np.empty(2)
                Hcit = h * cit3 / k4
                H2cit = h * Hcit / k3
                H3cit = H2cit * h / k2
                ni2pfree = (nio2 * k1 * (h ** (2))) / (1 + ((nio2 * k1 * (h ** (2))) / ni_total))
                NiH2cit = k7 * ni2pfree * H2cit
                NiHcit = k6 * ni2pfree * Hcit
                Nicit = k5 * ni2pfree * cit3
                F[0] = citrate_total - Hcit - H2cit - H3cit - Nicit - NiHcit - NiH2cit - cit3
                F[1] = ni_total - Nicit - NiHcit - NiH2cit - ni2pfree - nio2
                return F
            z = fsolve(f, [0.1, 0.1])
            cit3 = z[0]
            nio2 = z[1]
            ni2pfree = (nio2 * k1 * (h ** (2))) / (1 + ((nio2 * k1 * (h ** (2))) / ni_total))
            Hcit = h * cit3 / k4
            H2cit = h * Hcit / k3
            H3cit = H2cit * h / k2
            NiH2cit = k7 * ni2pfree * H2cit
            NiHcit = k6 * ni2pfree * Hcit
            Nicit = k5 * ni2pfree * cit3
        elif citrate_total == 0:
            f = 10 ** (logk-2*pH_x)
            nip2free = f / (1 + f / nitotal)
            nio2 = nitotal - nip2free
            cit3 = Hcit =H2cit =H3cit =NiHcit =Nicit =NiH2cit=0
        return [cit3, nio2, ni2pfree, Hcit, H2cit, H3cit, NiH2cit, NiHcit, Nicit]

    cit3freeplot = []
    nio2freeplot = []
    ni2pfreeplot = []
    Hcitfreeplot = []
    H2citfreeplot = []
    H3citfreeplot = []
    NiH2citfreeplot = []
    NiHcitfreeplot = []
    Nicitfreeplot = []

    for pHval in pH_x:
        cit3freeplot.append(concs(citrate_total, ni_total, pHval)[0])
        nio2freeplot.append(concs(citrate_total, ni_total, pHval)[1])
        ni2pfreeplot.append(concs(citrate_total, ni_total, pHval)[2])
        Hcitfreeplot.append(concs(citrate_total, ni_total, pHval)[3])
        H2citfreeplot.append(concs(citrate_total, ni_total, pHval)[4])
        H3citfreeplot.append(concs(citrate_total, ni_total, pHval)[5])
        NiH2citfreeplot.append(concs(citrate_total, ni_total, pHval)[6])
        NiHcitfreeplot.append(concs(citrate_total, ni_total, pHval)[7])
        Nicitfreeplot.append(concs(citrate_total, ni_total, pHval)[8])
    #----------------------------------------------------------------------------------------------
    # previous section generates the concentrations of species, in ni(oh)2 region complexes need to
    # only show if they are greater than ni(oh)2 conc, (n is for nh3, nhp for nh4+),
    # note pH_x input needs to be a list also
    NiH2cit = NiHcit = Nicit = nip2 =ni_total
    cit3 = Hcit = H2cit = H3cit =citrate_total
    T_ = 298

    def trace_generator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_):
        def interceptgenerator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_):
            if citrate_total != 0:
                bottom_eqs = [[R1(nip2, T_) for i in range(len(pH_x))], R2(pH_x, T_, NiH2cit, H3cit),
                              R3(pH_x, T_, NiHcit, H3cit),
                              R4(pH_x, T_, NiHcit, H2cit), R5(pH_x, T_, Nicit, H2cit), R6(pH_x, T_, Nicit, Hcit),
                              [R7(T_, Nicit, cit3) for i in range(len(pH_x))], R8(pH_x, T_)]
            elif citrate_total == 0:
                bottom_eqs = [[R1(nip2, T_) for i in range(len(pH_x))], R8(pH_x, T_)]

            y_vars = []
            # y vars has m, +c of every line in bottom eqs
            for i, eq in enumerate(bottom_eqs):
                y_vars.append(np.polyfit(pH_x, eq, 1))  # 1 means linear equaiton

            As = []  # M
            Bs = []  # C
            for i, var in enumerate(y_vars):
                if i >= 1:
                    As.append(np.array([[-y_vars[i - 1][0], 1], [-y_vars[i][0], 1]]))
                    Bs.append(np.array([y_vars[i - 1][1], y_vars[i][1]]))
            # inters has [x_intercept, y_intercept] of all intercepts
            inters = []
            for i, ms in enumerate(As):
                for j, cs in enumerate(Bs):
                    if i == j:
                        # inters.append(np.linalg.solve(As[i],Bs[j]))
                        inters.append(np.linalg.inv(As[i]).dot(Bs[j]))
            return inters

        inters = interceptgenerator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)

        xinters = []
        for item in inters:
            xinters.append(item[0])

        x_data = []
        if citrate_total != 0:
            for i, item in enumerate(xinters):
                if i == 0:
                    x_data.append(np.linspace(0, item, 5))
                elif i >= 1:
                    x_data.append(np.linspace(xinters[i - 1], item, 5))
            finalindex = len(xinters) - 1
            x_data.append(np.linspace(xinters[finalindex], 14, 5))
        elif citrate_total == 0:
            x_data.append(list(np.linspace(0, xinters[0], 5)))
            x_data.append(list(np.linspace(xinters[0], 14, 5)))

        if citrate_total != 0:
            y_data_bottom = [[R1(nip2, T_) for i in range(len(x_data[0]))], R2(x_data[1], T_, NiH2cit, H3cit),
                             R3(x_data[2], T_, NiHcit, H3cit),
                             R4(x_data[3], T_, NiHcit, H2cit), R5(x_data[4], T_, Nicit, H2cit),
                             R6(x_data[5], T_, Nicit, Hcit),
                             [R7(T_, Nicit, cit3) for i in range(len(x_data[6]))], R8(x_data[7], T_)]
        if citrate_total == 0:
            y_data_bottom = [[R1(nip2, T_) for i in range(len(x_data[0]))], R8(x_data[1], T_)]
        new_x_bottom = []
        new_y_bottom = []

        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_bottom.append(xvalue)
        for yvalues in y_data_bottom:
            for yvalue in yvalues:
                new_y_bottom.append(yvalue)

        if citrate_total != 0:
            y_data_top = [T1(nip2, x_data[0], T_), T2(x_data[1], NiH2cit, H3cit, T_), T3(x_data[2], NiHcit, H3cit, T_),
                          T4(x_data[3], NiHcit, H2cit, T_),
                          T5(x_data[4], Nicit, H2cit, T_), T6(x_data[5], Nicit, Hcit, T_),
                          T7(x_data[6], Nicit, cit3, T_), T8(x_data[7], T_)]
        elif citrate_total == 0:
            y_data_top = [T1(nip2, x_data[0], T_), T8(x_data[1], T_)]
        new_x_top = []
        new_y_top = []
        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_top.append(xvalue)
        for yvalues in y_data_top:
            for yvalue in yvalues:
                new_y_top.append(yvalue)

        if citrate_total != 0:
            y_interps = [T1(nip2, inters[0][0], T_), T2(inters[1][0], NiH2cit, H3cit, T_),
                         T3(inters[2][0], NiHcit, H3cit, T_), T4(inters[3][0], NiHcit, H2cit, T_),
                         T5(inters[4][0], Nicit, H2cit, T_), T6(inters[5][0], Nicit, Hcit, T_),
                         T7(inters[6][0], Nicit, cit3, T_)]
        elif citrate_total == 0:
            y_interps = [T1(nip2, inters[0][0], T_)]
        vys = []
        for i, val in enumerate(inters):
            vys.append(list(np.linspace(inters[i][1], y_interps[i], 5)))
        x_data_verticals = []
        for i in range(len(inters)):
            x_data_verticals.append([inters[i][0] for j in range(len(vys[i]))])
        new_x_vert = []
        new_y_vert = []
        for xvalues in x_data_verticals:
            for xvalue in xvalues:
                new_x_vert.append(xvalue)
        for yvalues in vys:
            for yvalue in yvalues:
                new_y_vert.append(yvalue)

        nio3regionx = list(new_x_bottom) + list([14 for i in range(0, 5)]) + list(
            reversed(np.linspace(0, 14, 5))) + list(
            [0 for i in range(0, 5)])
        nio3regiony = list(new_y_top) + list(np.linspace(T8(14, T_), 2.6, 5)) + list([2.6 for i in range(0, 5)]) + list(
            np.linspace(T1(nip2, 0, 298), 2.6, 5))

        niregionx = list(new_x_bottom) + list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list(
            [0 for i in range(0, 5)])
        niregiony = list(new_y_bottom) + list(np.linspace(R8(14, T_), -1.8, 5)) + list(
            [-1.8 for i in range(0, 5)]) + list(
            np.linspace(-1.8, R1(nip2, 298), 5))

        nip2regionx = list(x_data[0]) + list(x_data_verticals[0]) + list(reversed(x_data[0])) + list(
            [0 for i in range(0, 5)])
        nip2regiony = list(y_data_bottom[0]) + list(vys[0]) + list(reversed(y_data_top[0])) + list(
            np.linspace(R1(nip2, T_), T1(nip2, 0, 298), 5))

        if citrate_total != 0:
            NiH2citregionx = list(x_data[1]) + list(x_data[2]) + list(x_data_verticals[0]) + list(x_data[1]) + list(
                x_data[2]) + list(x_data_verticals[1])
            NiH2citregiony = list(y_data_bottom[1]) + list(y_data_bottom[2]) + list(vys[0]) + list(
                y_data_top[1]) + list(y_data_top[2]) + list(reversed(vys[1]))

            NiHcitregionx = list(x_data[3]) + list(x_data_verticals[1]) + list(x_data[3]) + list(x_data_verticals[3])
            NiHcitregiony = list(y_data_bottom[3]) + list(vys[1]) + list(y_data_top[3]) + list(reversed(vys[3]))

            Nicitregionx = list(x_data[4]) + list(x_data[5]) + list(x_data[6]) + list(x_data_verticals[3]) + list(
                x_data[4]) + list(x_data[5]) + list(x_data[6]) + list(x_data_verticals[6])
            Nicitregiony = list(y_data_bottom[4]) + list(y_data_bottom[5]) + list(y_data_bottom[6]) + list(
                vys[3]) + list(y_data_top[4]) + list(y_data_top[5]) + list(y_data_top[6]) + list(reversed(vys[6]))

            nio2regionx = list(reversed(x_data[7])) + list(x_data_verticals[6]) + list(x_data[7]) + list(
                [14 for i in range(0, 5)])
            nio2regiony = list(reversed(y_data_bottom[7])) + list(vys[6]) + list(y_data_top[7]) + list(
                np.linspace(R8(14, T_), T8(14, T_), 5))
            xs = [nip2regionx, NiH2citregionx, NiHcitregionx, Nicitregionx, nio2regionx]
            ys = [nip2regiony, NiH2citregiony, NiHcitregiony, Nicitregiony, nio2regiony]

        elif citrate_total == 0:
            nio2regionx = list(reversed(x_data[1])) + list(x_data_verticals[0]) + list(x_data[1])
            nio2regiony = list(reversed(y_data_bottom[1])) + list(vys[0]) + list(y_data_top[1])
            xs = [nip2regionx, nio2regionx]
            ys = [nip2regiony, nio2regiony]

        return [xs, ys, niregionx, niregiony, nio3regionx, nio3regiony]

    # end function, return data to add to traces, should already be in correct form to
    # cooperate with dash notation

    xs = trace_generator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[0]
    ys = trace_generator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[1]
    niregionx = trace_generator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[2]
    niregiony = trace_generator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[3]
    nio3regionx = trace_generator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[4]
    nio3regiony = trace_generator(pH_x, nip2, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[5]

    if citrate_total != 0:
        name= ['Ni<sup>2+</sup>', 'NiH<sub>2</sub>cit<sup>+</sup>',  'NiHcit</sub>', 'Nicit<sup>-</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)',  'rgba(252, 177, 101, 0.8)', 'rgba(7, 117, 189, 0.66)', 'rgba(63, 63, 191, 0.5)']
    elif citrate_total == 0:
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        color= ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)']

    data = []
    if citrate_total != 0:
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

    elif citrate_total == 0:
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
    extrax = [niregionx, nio3regionx]
    extray = [niregiony, nio3regiony]
    namee = ['Ni', 'Ni(OH)<sub>3</sub>']
    colors1 = ['rgba(127, 63, 191, 0.5)', 'rgba(30, 205, 40, 0.5)']
    tracesextra = []
    for i, extraxx in enumerate(extrax):
        tracesextra.append(go.Scatter(
            x=extraxx,
            y=extray[i],
            mode='lines',
            name=namee[i],
            fill='toself',
            showlegend=True,
            hoverinfo='skip'
        ))

    for trace in tracesextra:
        fig.add_traces(trace)
    if citrate_total != 0:
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
            x=2.1,
            y=2.2,
            text='H<sub>3</sub>cit',
            showarrow=False,
            font=dict(
                    family="Courier New bold",
                    size=16,
                    color="purple"
                    ),)
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
####################################################################################################################
@app.callback(
    Output('speciation plot', 'figure'),
    [Input('nickel_slider', 'value'),
     Input('citrate_dropdown', 'value')])

def speciation_graph(citrate_total, ni_total):
    k1 = 9.197 * (10 ** 11)
    k2 = 1.01 * (10 ** -3)
    k3 = 1.787 * (10 ** -5)
    k4 = 4.031 * (10 ** -7)
    k5 = 2.512 * (10 ** 5)
    k6 = 1.995 * (10 ** 3)
    k7 = 5.623 * (10 ** 1)
    pH_x = np.linspace(0, 14, 71)
    #----------------------------------------------------------------------------------------------
    # begin first function, output all species concentrations. One concentration for each pH value.

    def concs(citrate_total, ni_total, pH_x):
        h = 10 ** (-pH_x)

        def f(z):
            cit3 = z[0]
            nio2 = z[1]
            F = np.empty(2)
            Hcit = h * cit3 / k4
            H2cit = h * Hcit / k3
            H3cit = H2cit * h / k2
            ni2pfree = (nio2 * k1 * (h ** (2))) / (1 + ((nio2 * k1 * (h ** (2))) / ni_total))
            NiH2cit = k7 * ni2pfree * H2cit
            NiHcit = k6 * ni2pfree * Hcit
            Nicit = k5 * ni2pfree * cit3
            F[0] = citrate_total - Hcit - H2cit - H3cit - Nicit - NiHcit - NiH2cit - cit3
            F[1] = ni_total - Nicit - NiHcit - NiH2cit - ni2pfree - nio2
            return F
        z = fsolve(f, [0.1, 0.1])
        cit3 = z[0]
        nio2 = z[1]
        ni2pfree = (nio2 * k1 * (h ** (2))) / (1 + ((nio2 * k1 * (h ** (2))) / ni_total))
        Hcit = h * cit3 / k4
        H2cit = h * Hcit / k3
        H3cit = H2cit * h / k2
        NiH2cit = k7 * ni2pfree * H2cit
        NiHcit = k6 * ni2pfree * Hcit
        Nicit = k5 * ni2pfree * cit3
        # nio2_1=ni2pfree/(k1*(h**2))
        # ni2pfree_1= ni_total-Nicit-NiHcit-NiH2cit-nio2_1
        return [cit3, nio2, ni2pfree, Hcit, H2cit, H3cit, NiH2cit, NiHcit, Nicit]


    cit3freeplot = []
    nio2freeplot = []
    ni2pfreeplot = []
    Hcitfreeplot = []
    H2citfreeplot = []
    H3citfreeplot = []
    NiH2citfreeplot = []
    NiHcitfreeplot = []
    Nicitfreeplot = []

    for pHval in pH_x:
        cit3freeplot.append(concs(citrate_total, ni_total, pHval)[0])
        nio2freeplot.append(concs(citrate_total, ni_total, pHval)[1])
        ni2pfreeplot.append(concs(citrate_total, ni_total, pHval)[2])
        Hcitfreeplot.append(concs(citrate_total, ni_total, pHval)[3])
        H2citfreeplot.append(concs(citrate_total, ni_total, pHval)[4])
        H3citfreeplot.append(concs(citrate_total, ni_total, pHval)[5])
        NiH2citfreeplot.append(concs(citrate_total, ni_total, pHval)[6])
        NiHcitfreeplot.append(concs(citrate_total, ni_total, pHval)[7])
        Nicitfreeplot.append(concs(citrate_total, ni_total, pHval)[8])

    datasets = [cit3freeplot, nio2freeplot, ni2pfreeplot, Hcitfreeplot, H2citfreeplot,
                H3citfreeplot, NiH2citfreeplot, NiHcitfreeplot, Nicitfreeplot]
    name = ['Cit<sup>3-</sup>', 'Ni(OH)<sub>2</sub>', 'Ni<sup>2+</sup>', 'Hcit<sup>2-</sup>',
            'H<sub>2</sub>cit<sup>2-</sup>','H<sub>3</sub>cit', 'NiH<sub>2</sub>cit<sup>+</sup>',
            'NiH<sub>cit</sub>', 'Nicit<sup>-</sup>']
    fill = [None, None, None, None, None, None, None, None, None]
    color = ['rgb(90, 0, 100)', 'rgb(40, 130, 80)', 'rgb(9, 0, 0)', 'rgb(63, 63, 191)', 'rgb(191, 63, 63)',
             'rgb(66, 81, 245)', 'rgb(218, 66, 245)', 'rgb(245, 144, 66)', 'rgb(245, 66, 90)']
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
                    range=[0, max(citrate_total, ni_total)*1.05])
    fig1.update_layout(
        title={
            'text': "Speciation plot",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig1
