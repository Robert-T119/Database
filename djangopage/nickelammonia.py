from djangopage.reactions import *
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import fsolve
from django_plotly_dash import DjangoDash

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app = DjangoDash('nickelammonia', add_bootstrap_links=True)
# , add_bootstrap_links=True
# add_bootstrap_links=True


app.layout = html.Div([
    html.Div([
        html.H2(u'Nickel-Ammonia-H\u2082O system')],
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
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]},
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
        html.H6(u'Total ammonia concentration (kmolm\u207B\u00B3):'),
        dcc.Slider(
            id='ammonia_dropdown',
            min=0.0,
            max=3.0,
            value=1.2,
            step=0,
            marks={i: str(i) for i in [0, 0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]},
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
     Input('ammonia_dropdown', 'value')])

def speciation_graph(nitotal, ntotal):
    T_ = 298
    k2 = 10 ** (-9.25)
    k4 = 10 ** (-18.26)
    k5 = 10 ** (-35.91)
    logk = 11.96
    K0 = 5.65 * 10 ** (-10)
    #----------------------------------------------------------------------------------------------
    # begin first function, output all species concentrations. One concentration for each pH value.
    def species1(pH_x, ntotal, nitotal):
        h = 10**(-pH_x)
        if ntotal != 0.0:
            f = 10**(logk - 2 * pH_x)
            nip2free = f/(1+f/nitotal)
            h = 10**(-pH_x)
            def solve1(nhp):
                nin4p2 = k4*(nhp**(4))/((h)**(2))
                nin6p2 = k5*(nhp**(6))/((h)**(4))
                return (ntotal-nhp-k2*nhp/h-4*nin4p2-6*nin6p2)
            nhp = fsolve(solve1,0.5)
            nh3 = k2*nhp/h
            nin4p2 =k4*(nhp**(4))/((h)**(2))
            nin6p2 = k5*(nhp**(6))/((h)**(4))
            nio2 = nitotal-nin4p2-nin6p2-nip2free
        elif ntotal == 0.0:
            f = 10**(logk - 2 * pH_x)
            nip2free = f/(1+f/nitotal)
            nio2 = nitotal - nip2free
            nin4p2 = nin6p2 = nh3 = nhp = 0
        return [nh3, nhp, nin4p2, nin6p2, nip2free, nio2]

    pH_x = list(np.linspace(0, 14, 50))
    nplot = []
    nhpplot = []
    nin4p2plot = []
    nin6p2plot = []
    nip2plot = []
    nip2pptplot = []
    for pHval in pH_x:
        nplot.append(species1(pHval, ntotal, nitotal)[0])
        nhpplot.append(species1(pHval, ntotal, nitotal)[1])
        nin4p2plot.append(species1(pHval, ntotal, nitotal)[2])
        nin6p2plot.append(species1(pHval, ntotal, nitotal)[3])
        nip2plot.append(species1(pHval, ntotal, nitotal)[4])
        nip2pptplot.append(species1(pHval, ntotal, nitotal)[5])
    #----------------------------------------------------------------------------------------------
    # previous section generates the concentrations of species, in ni(oh)2 region complexes need to
    # only show if they are greater than ni(oh)2 conc, (n is for nh3, nhp for nh4+),
    # note pH_x input needs to be a list also
    nip2 = nip2ppt = nitotal
    nhp = n = ntotal
    if nip2 <= nhp / 4:
        nin4p2 = nip2
    else:
        nin4p2 = nhp / 4
    if nip2 <= nhp / 6:
        nin6p2 = nip2
    else:
        nin6p2 = nhp / 6

    def nio2checker(nin4p2plot, nin6p2plot, nip2pptplot, n):
        if n != 0.0:
            maxnin4p2 = max(nin4p2plot)
            i = nin4p2plot.index(maxnin4p2)
            maxnin6p2 = max(nin6p2plot)
            j = nin6p2plot.index(maxnin6p2)
            if maxnin6p2 > nip2pptplot[j] and maxnin4p2 > nip2pptplot[i]:
                status = 0
            elif maxnin6p2 > nip2pptplot[j] and maxnin4p2 < nip2pptplot[i]:
                status = 1
            elif maxnin6p2 < nip2pptplot[j] and maxnin4p2 > nip2pptplot[i]:
                status = 2
            elif maxnin6p2 < nip2pptplot[j] and maxnin4p2 < nip2pptplot[i]:
                status = 3
        elif n == 0.0:
            status = 3
        return status

    status = nio2checker(nin4p2plot, nin6p2plot, nip2pptplot, n)
    # nio2checker(nin4p2plot, nin6p2plot, nip2pptplot)
    # nio2checker returns status. This can take 4 values, depending on what will be shown inside
    # the ni(oh)2 region
    # --------------------------------------------------------------------------------------------------------------
    def trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_):
        def interceptgenerator(pH_x, nip2, nin4p2, nin6p2, n, nhp, T_):
            if status == 0:
                bottom_eqs = [[R1(nip2, T_) for i in range(len(pH_x))], R2(pH_x, T_), R3(pH_x, nhp, nin4p2, T_),
                      R4(pH_x, nhp, nin6p2, T_), [R6(nhp, nin6p2, T_) for i in range(len(pH_x))]]
            elif status == 1:
                bottom_eqs = [[R1(nip2, T_) for i in range(len(pH_x))], R2(pH_x, T_),
                              R4(pH_x, nhp, nin6p2, T_), [R6(nhp, nin6p2, T_) for i in range(len(pH_x))]]
            elif status == 2:
                bottom_eqs = [[R1(nip2, T_) for i in range(len(pH_x))], R2(pH_x, T_), R3(pH_x, nhp, nin4p2, T_)]
            elif status == 3 or n == 0.0:
                bottom_eqs = [[R1(nip2, T_) for i in range(len(pH_x))], R2(pH_x, T_)]

            y_vars = []
            # y vars has m, +c of every line in bottom eqs
            for i, eq in enumerate(bottom_eqs):
                y_vars.append(np.polyfit(pH_x, eq, 1))
            if status != 3:
                y_vars.append(np.polyfit(pH_x, R2(pH_x, T_), 1))
            As = []
            Bs = []
            for i, var in enumerate(y_vars):
                if i >= 1:
                    As.append(np.array([[-y_vars[i-1][0], 1], [-y_vars[i][0], 1]]))
                    Bs.append(np.array([y_vars[i-1][1], y_vars[i][1]]))
            # inters has [x_intercept, y_intercept] of all intercepts
            inters = []
            for i, ms in enumerate(As):
                for j, cs in enumerate(Bs):
                    if i == j:
                        inters.append(np.linalg.inv(As[i]).dot(Bs[j]))
            return inters
        # interseptgenerator returns inters is array of x and y values of intercepts
        #------------------------------------------------------------------------------------
        # inters = interceptgenerator(pH_x, nip2, nin4p2, nin6p2, n, nhp, T_, status)[0]
        inters = interceptgenerator(pH_x, nip2, nin4p2, nin6p2, n, nhp, T_)
        xinters = []
        for item in inters:
            xinters.append(item[0])
        x_data = []
        if status != 3:
            for i, item in enumerate(xinters):
                if i == 0:
                    x_data.append(np.linspace(0, item, 5))
                elif i >= 1:
                    x_data.append(np.linspace(xinters[i-1], item, 5))
            finalindex = len(xinters)-1
            x_data.append(np.linspace(xinters[finalindex], 14, 5))
        elif status == 3 or n == 0.0:
            x_data.append(list(np.linspace(0, xinters[0], 5)))
            x_data.append(list(np.linspace(xinters[0], 14, 5)))

        if status == 0:
            y_data_bottom = [[R1(nip2, T_) for i in range(len(x_data[0]))], R2(x_data[1], T_), R3(x_data[2], nhp, nin4p2, T_),
                             R4(x_data[3], nhp, nin6p2, T_), [R6(nhp, nin6p2, T_) for i in range(len(x_data[4]))],
                             R2(x_data[5], T_)]
        elif status == 1:
            y_data_bottom = [[R1(nip2, T_) for i in range(len(x_data[0]))], R2(x_data[1], T_), R4(x_data[2], nhp, nin6p2, T_),
                             [R6(nhp, nin6p2, T_) for i in range(len(x_data[3]))], R2(x_data[4], T_)]
        elif status == 2:
            y_data_bottom = [[R1(nip2, T_) for i in range(len(x_data[0]))], R2(x_data[1], T_), R3(x_data[2], nhp, nin4p2, T_),
                             R2(x_data[3], T_)]
        elif status == 3 or n == 0.0:
            y_data_bottom = [[R1(nip2, T_) for i in range(len(x_data[0]))], R2(x_data[1], T_)]

        new_x_bottom = []
        new_y_bottom = []
        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_bottom.append(xvalue)
        for yvalues in y_data_bottom:
            for yvalue in yvalues:
                new_y_bottom.append(yvalue)

        if status == 0:
            y_data_top = [T1(nip2, x_data[0], T_), T2(x_data[1], T_), T3(nhp, nin4p2, x_data[2], T_),
                          T4(x_data[3], nhp, nin6p2, T_), T5(nhp, x_data[4], nin6p2, T_), T2(x_data[5], T_)]
        elif status == 1:
            y_data_top = [T1(nip2, x_data[0], T_), T2(x_data[1], T_), T4(x_data[2], nhp, nin6p2, T_),
                          T5(nhp, x_data[3], nin6p2, T_), T2(x_data[4], T_)]
        elif status == 2:
            y_data_top = [T1(nip2, x_data[0], T_), T2(x_data[1], T_), T3(nhp, nin4p2, x_data[2], T_),
                          T2(x_data[3], T_)]
        elif status == 3 or n == 0.0:
            y_data_top = [T1(nip2, x_data[0], T_), T2(x_data[1], T_)]
        new_x_top = []
        new_y_top = []
        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_top.append(xvalue)
        for yvalues in y_data_top:
            for yvalue in yvalues:
                new_y_top.append(yvalue)

        if status == 0:
            y_interps = [T1(nip2, inters[0][0], T_), T2(inters[1][0], T_), T3(nhp, nin4p2, inters[2][0], T_),
                         T4(inters[3][0], nhp, nin6p2, T_), T5(nhp, inters[4][0], nin6p2, T_)]
        elif status == 1:
            y_interps = [T1(nip2, inters[0][0], T_), T2(inters[1][0], T_), T4(inters[2][0], nhp, nin6p2, T_),
                         T5(nhp, inters[3][0], nin6p2, T_)]
        elif status == 2:
            y_interps = [T1(nip2, inters[0][0], T_), T2(inters[1][0], T_), T3(nhp, nin4p2, inters[2][0], T_)]
        elif status == 3 or n == 0.0:
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

        nip2regionx = list(x_data[0]) + list(x_data_verticals[0]) + list(reversed(x_data[0]))
        nip2regiony = list(y_data_bottom[0]) + list(vys[0]) + list(reversed(y_data_top[0]))
        nio3regiony = list(new_y_top) + list(np.linspace(T2(14, T_), 2.4, 5)) + list([2.4 for i in range(0, 5)]) + list(np.linspace(T1(nitotal, 0, 298), 2.4, 5))
        nio3regionx = list(new_x_bottom) + list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list([0 for i in range(0, 5)])
        niregiony = list(new_y_bottom) + list(np.linspace(R2(14, T_), -1, 5)) + list([-1 for i in range(0, 5)]) + list(np.linspace(-1, R1(nitotal, 298), 5))
        niregionx = list(new_x_bottom) + list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list([0 for i in range(0, 5)])
        if status == 0:
            nio2regionx = list(x_data[1]) + list(x_data_verticals[0]) + list(x_data[1]) + list(x_data_verticals[1])
            nio2regiony = list(y_data_bottom[1])  + list(vys[0]) + list(y_data_top[1]) + list(reversed(vys[1]))
            nin4p2regionx = list(x_data[2]) + list(x_data_verticals[1]) + list(x_data[2]) + list(x_data_verticals[2])
            nin4p2regiony = list(y_data_bottom[2]) + list(vys[1]) + list(y_data_top[2]) + list(reversed(vys[2]))
            nin6p2regionx = list(reversed(x_data[3])) + list(reversed(x_data[4])) + list(x_data_verticals[2]) + list(x_data[3]) + list(x_data[4]) + list(x_data_verticals[4])
            nin6p2regiony = list(reversed(y_data_bottom[3])) + list(reversed(y_data_bottom[4])) + list(vys[2]) + list(y_data_top[3]) + list(y_data_top[4]) + list(reversed(vys[4]))
            unlabelledx = list(reversed(x_data[5])) + list(x_data_verticals[4]) + list(x_data[5])
            unlabelledy = list(reversed(y_data_bottom[5])) + list(vys[4]) + list(y_data_top[5])
            data = [nip2regionx, nip2regiony, nio2regionx, nio2regiony, unlabelledx, unlabelledy, nin4p2regionx, nin4p2regiony,
               nin6p2regionx, nin6p2regiony]
            xs = [data[i] for i in np.arange(0, 10, 2)]
            ys = [data[i] for i in np.arange(1, 11, 2)]
        elif status == 1:
            nio2regionx = list(x_data[1]) + list(x_data_verticals[0]) + list(x_data[1]) + list(x_data_verticals[1])
            nio2regiony = list(y_data_bottom[1])  + list(vys[0]) + list(y_data_top[1]) + list(reversed(vys[1]))
            nin6p2regionx = list(reversed(x_data[2])) + list(reversed(x_data[3])) + list(x_data_verticals[1]) + list(x_data[2]) + list(x_data[3]) + list(x_data_verticals[3])
            nin6p2regiony = list(reversed(y_data_bottom[2])) + list(reversed(y_data_bottom[3])) + list(vys[1]) + list(y_data_top[2]) + list(y_data_top[3]) + list(reversed(vys[3]))
            unlabelledx = list(reversed(x_data[4])) + list(x_data_verticals[3]) + list(x_data[4])
            unlabelledy = list(reversed(y_data_bottom[4])) + list(vys[3]) + list(y_data_top[4])
            xs = [nip2regionx, nio2regionx, unlabelledx, nin6p2regionx]
            ys = [nip2regiony, nio2regiony, unlabelledy, nin6p2regiony]
        elif status == 2:
            nio2regionx = list(x_data[1]) + list(x_data_verticals[0]) + list(x_data[1]) + list(x_data_verticals[1])
            nio2regiony = list(y_data_bottom[1])  + list(vys[0]) + list(y_data_top[1]) + list(reversed(vys[1]))
            nin4p2regionx = list(x_data[2]) + list(x_data_verticals[1]) + list(x_data[2]) + list(x_data_verticals[2])
            nin4p2regiony = list(y_data_bottom[2]) + list(vys[1]) + list(y_data_top[2]) + list(reversed(vys[2]))
            unlabelledx = list(reversed(x_data[3])) + list(x_data_verticals[2]) + list(x_data[3])
            unlabelledy = list(reversed(y_data_bottom[3])) + list(vys[2]) + list(y_data_top[3])
            xs = [nip2regionx, nio2regionx, unlabelledx, nin4p2regionx]
            ys = [nip2regiony, nio2regiony, unlabelledy, nin4p2regiony]
        elif status == 3:
            nio2regionx =  list(reversed(x_data[1])) + list(x_data_verticals[0]) + list(x_data[1])
            nio2regiony = list(reversed(y_data_bottom[1]))  + list(vys[0]) + list(y_data_top[1])
            xs = [nip2regionx, nio2regionx]
            ys = [nip2regiony, nio2regiony]
        return [xs, ys, niregionx, niregiony, nio3regionx, nio3regiony]
    # end function, return data to add to traces, should already be in correct form to
    # cooperate with dash notation
    #------------------------------------------------------------------------------------------------
    xs = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[0]
    ys = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[1]
    niregionx = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[2]
    niregiony = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[3]
    nio3regionx = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[4]
    nio3regiony = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[5]
    if status == 0:
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>',  'Ni(OH)<sub>2</sub>', '[Ni(NH<sub>3</sub>)<sub>4</sub>]<sup>2+</sup>', '[Ni(NH<sub>3</sub>)<sub>6</sub>]<sup>2+</sup>']
        color = ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)',  'rgba(191, 191, 63, 0.5)', 'rgba(9, 0, 0, 0.9)', 'rgba(63, 63, 191, 0.5)']
    elif status == 1:
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>', 'Ni(OH)<sub>2</sub>', '[Ni(NH<sub>3</sub>)<sub>6</sub>]<sup>2+</sup>']
        color = ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)',  'rgba(191, 191, 63, 0.5)', 'rgba(63, 63, 191, 0.5)']
    elif status == 2: # this never really happens
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>', 'Ni(OH)<sub>2</sub>', '[Ni(NH<sub>3</sub>)<sub>4</sub>]<sup>2+</sup>']
    elif status == 3:
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)']
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
    extrax = [niregionx, nio3regionx]
    extray = [niregiony, nio3regiony]
    namee = ['Ni', 'Ni(OH)<sub>3</sub>']
    colors1 = ['rgba(127, 63, 191, 0.5)', 'rgba(30, 205, 40, 0.5)']
    tracesextra = []
    for i, extraxx in enumerate(extrax):
        tracesextra.append(go.Scatter(
            x = extraxx,
            y = extray[i],
            mode='lines',
            name=namee[i],
            fill='toself',
            showlegend=True,
            hoverinfo='skip'
        ))
    for trace in tracesextra:
        fig.add_traces(trace)
    if ntotal != 0:
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
    Input('ammonia_dropdown', 'value')])

def speciation_graph(nitotal, ntotal):
    k2 = 10 ** (-9.25)
    k4 = 10 ** (-18.26)
    k5 = 10 ** (-35.91)
    logk = 11.96
    K0 = 5.65 * 10 ** (-10)
    #----------------------------------------------------------------------------------------------
    # begin first function, output all species concentrations. One concentration for each pH value.
    def species1(pH_x, ntotal, nitotal):
        h = 10**(-pH_x)
        if ntotal != 0.0:
            f = 10**(logk - 2 * pH_x)
            nip2free = f/(1+f/nitotal)
            h = 10**(-pH_x)
            def solve1(nhp):
                nin4p2 = k4*(nhp**(4))/((h)**(2))
                nin6p2 = k5*(nhp**(6))/((h)**(4))
                return (ntotal-nhp-k2*nhp/h-4*nin4p2-6*nin6p2)
            nhp = fsolve(solve1,0.5)
            nh3 = k2*nhp/h
            nin4p2 =k4*(nhp**(4))/((h)**(2))
            nin6p2 = k5*(nhp**(6))/((h)**(4))
            nio2 = nitotal-nin4p2-nin6p2-nip2free
        elif ntotal == 0.0:
            f = 10**(logk - 2 * pH_x)
            nip2free = f/(1+f/nitotal)
            nio2 = nitotal - nip2free
            nin4p2 = nin6p2 = nh3 = nhp = 0
        return [nh3, nhp, nin4p2, nin6p2, nip2free, nio2]

    pH_x = list(np.linspace(0, 14, 50))
    nplot = []
    nhpplot = []
    nin4p2plot = []
    nin6p2plot = []
    nip2plot = []
    nip2pptplot = []
    for pHval in pH_x:
        nplot.append(species1(pHval, ntotal, nitotal)[0])
        nhpplot.append(species1(pHval, ntotal, nitotal)[1])
        nin4p2plot.append(species1(pHval, ntotal, nitotal)[2])
        nin6p2plot.append(species1(pHval, ntotal, nitotal)[3])
        nip2plot.append(species1(pHval, ntotal, nitotal)[4])
        nip2pptplot.append(species1(pHval, ntotal, nitotal)[5])
    if ntotal != 0.0:
        datasets = [[i[0] for i in nplot], [i[0] for i in nhpplot], [i[0] for i in nin4p2plot],
                [i[0] for i in nin6p2plot], nip2plot, [i[0] for i in nip2pptplot]]
        name = ['NH<sub>3</sub>', 'NH<sub>4</sub><sup>+</sup>', '[Ni(NH<sub>3</sub>)<sub>4</sub>]<sup>2+</sup>', '[Ni(NH<sub>3</sub>)<sub>6</sub>]<sup>2+</sup>', 'Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        fill = [None, None, None, None, None, None]
        color = ['rgb(90, 0, 100)', 'rgb(40, 130, 80)', 'rgb(245, 137, 22)', 'rgb(63, 63, 191)', 'rgb(191, 63, 63)', 'rgb(15, 15, 15)']
    elif ntotal == 0.0:
        datasets = [nip2plot, nip2pptplot]
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        fill = [None for i in range(len(name))]
        color = ['rgb(191, 63, 63)', 'rgb(243, 238, 77)']

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
                    range=[0, max(nitotal, ntotal)*1.05])
    fig1.update_layout(
        title={
            'text': "Speciation plot",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig1
