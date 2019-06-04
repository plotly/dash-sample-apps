# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_colorscales as dcs
import numpy as  np
import json
from textwrap import dedent as d

app = dash.Dash()

DEFAULT_COLORSCALE = [[0, 'rgb(12,51,131)'], [0.25, 'rgb(10,136,186)'],\
        [0.5, 'rgb(242,211,56)'], [0.75, 'rgb(242,143,56)'], \
        [1, 'rgb(217,30,30)']]

DEFAULT_COLORSCALE_NO_INDEX = [ea[1] for ea in DEFAULT_COLORSCALE]

def read_mniobj(file):
    ''' function to read a MNI obj file '''

    def triangulate_polygons(list_vertex_indices):
        ''' triangulate a list of n indices  n=len(list_vertex_indices) '''

        for k in range(0, len(list_vertex_indices), 3):
                yield list_vertex_indices[k: k+3]

    fp=open(file,'r')
    n_vert=[]
    n_poly=[]
    k=0
    list_indices=[]
    # Find number of vertices and number of polygons, stored in .obj file.
    # Then extract list of all vertices in polygons
    for i, line in enumerate(fp):
        if i==0:
        #Number of vertices
             n_vert=int(line.split()[6])
             vertices=np.zeros([n_vert,3])
        elif i<=n_vert:
             vertices[i-1]=list(map(float, line.split()))
        elif i>2*n_vert+5:
            if not line.strip():
                k=1
            elif k==1:
                list_indices.extend(line.split())
    #at this point list_indices is a list of strings, and each string is a vertex index, like this '23'
    #maps in Python 3.6 returns a generator, hence we convert it to a list
    list_indices=list(map(int, list_indices))#conver the list of string indices to int indices
    faces=np.array(list(triangulate_polygons(np.array(list_indices))))
    return vertices, faces

def standard_intensity(x,y,z):
    ''' color the mesh with a colorscale according to the values
                              of the vertices z-coordinates '''
    return z

def plotly_triangular_mesh(vertices, faces, intensities=None, colorscale="Viridis",
                           flatshading=False, showscale=False, reversescale=False, plot_edges=False):
    ''' vertices = a numpy array of shape (n_vertices, 3)
        faces = a numpy array of shape (n_faces, 3)
        intensities can be either a function of (x,y,z) or a list of values '''

    x,y,z=vertices.T
    I,J,K=faces.T

    if intensities is None:
        intensities = standard_intensity(x,y,z)

    if hasattr(intensities, '__call__'):
        intensity=intensities(x,y,z)#the intensities are computed  via a function,
                                    #that returns the list of vertices intensities
    elif  isinstance(intensities, (list, np.ndarray)):
        intensity=intensities#intensities are given in a list
    else:
        raise ValueError("intensities can be either a function or a list, np.array")

    mesh=dict(
        type='mesh3d',
        x=x, y=y, z=z,
        colorscale=colorscale,
        intensity= intensities,
        flatshading=flatshading,
        i=I, j=J, k=K,
        name='',
        showscale=showscale
    )

    mesh.update(lighting=dict( ambient= 0.18,
                                  diffuse= 1,
                                  fresnel=  0.1,
                                  specular= 1,
                                  roughness= 0.1,
                                  facenormalsepsilon=1e-6,
                                  vertexnormalsepsilon= 1e-12))

    mesh.update(lightposition=dict(x=100,
                                      y=200,
                                      z= 0))

    if  showscale is True:
            mesh.update(colorbar=dict(thickness=20, ticklen=4, len=0.75))

    if plot_edges is False: # the triangle sides are not plotted
        return  [mesh]
    else:#plot edges
        #define the lists Xe, Ye, Ze, of x, y, resp z coordinates of edge end points for each triangle
        #None separates data corresponding to two consecutive triangles
        tri_vertices= vertices[faces]
        Xe=[]
        Ye=[]
        Ze=[]
        for T in tri_vertices:
            Xe+=[T[k%3][0] for k in range(4)]+[ None]
            Ye+=[T[k%3][1] for k in range(4)]+[ None]
            Ze+=[T[k%3][2] for k in range(4)]+[ None]
        #define the lines to be plotted
        lines=dict(type='scatter3d',
                   x=Xe,
                   y=Ye,
                   z=Ze,
                   mode='lines',
                   name='',
                   line=dict(color= 'rgb(70,70,70)', width=1)
               )
        return [mesh, lines]

pts, tri=read_mniobj("surf_reg_model_both.obj")
intensities=np.loadtxt('aal_atlas.txt')
data=plotly_triangular_mesh(pts, tri, intensities,
                            colorscale=DEFAULT_COLORSCALE, flatshading=False,
                            showscale=False, reversescale=False, plot_edges=False)
data[0]['name'] = 'human_atlas'

axis_template = dict(
    showbackground=True,
    backgroundcolor="rgb(10, 10,10)",
    gridcolor="rgb(255, 255, 255)",
    zerolinecolor="rgb(255, 255, 255)")

plot_layout = dict(
         title = '',
         margin=dict(t=0,b=0,l=0,r=0),
         font=dict(size=12, color='white'),
         width=700,
         height=700,
         showlegend=False,
         plot_bgcolor='black',
         paper_bgcolor='black',
         scene=dict(xaxis=axis_template,
                    yaxis=axis_template,
                    zaxis=axis_template,
                    aspectratio=dict(x=1, y=1.2, z=1),
                    camera=dict(eye=dict(x=1.25, y=1.25, z=1.25)),
                    annotations=[]
                )
        )

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'padding': '10px',
        'marginBottom': '20px'
    },
    'graph': {
        'userSelect': 'none',
        'margin': 'auto'
    }
}

'''
~~~~~~~~~~~~~~~~
~~ APP LAYOUT ~~
~~~~~~~~~~~~~~~~
'''

app.layout = html.Div(children=[
    html.P(
        children=['''
            Click on the brain to add an annotation. \
            Drag the black corners of the graph to rotate. ''',
            html.A(
                children='GitHub',
                target= '_blank',
                href='https://github.com/plotly/dash-brain-surface-viewer',
                style={'color': '#F012BE'}
            ),
            '.'
        ]
    ),
    html.Div([
        html.P(
            children='Click colorscale to change:',
            style={'display':'inline-block', 'fontSize':'12px'}
        ),
        html.Div([
    	    dcs.DashColorscales(
    	        id='colorscale-picker',
                colorscale=DEFAULT_COLORSCALE_NO_INDEX
            )], style={'marginTop':'-15px', 'marginLeft':'-30px'}
	    ),
        html.Div([
            dcc.RadioItems(
                options=[
                    {'label': 'Cortical Thickness', 'value': 'human'},
                    {'label': 'Mouse Brain', 'value': 'mouse'},
                    {'label': 'Brain Atlas', 'value': 'human_atlas'},
                ],
                value='human_atlas',
                id='radio-options',
                labelStyle={'display': 'inline-block'}
            )
        ])
	]),
    dcc.Graph(
        id='brain-graph',
        figure={
            'data': data,
            'layout': plot_layout,
        },
        config={'editable': True, 'scrollZoom': False},
        style=styles['graph']
    ),
    html.Div([
        dcc.Markdown(d("""
            **Click Data**

            Click on points in the graph.
        """)),
        html.Pre(id='click-data', style=styles['pre']),
    ]),
    html.Div([
        dcc.Markdown(d("""
            **Relayout Data**

            Drag the graph corners to rotate it.
        """)),
        html.Pre(id='relayout-data', style=styles['pre']),
    ]),
    html.P(
        children=[
            'Dash/Python code on ',
            html.A(
                children='GitHub',
                target= '_blank',
                href='https://github.com/plotly/dash-brain-surface-viewer',
                style={'color': '#F012BE'}
            ),
            '. Brain data from Mcgill\'s ACE Lab  ',
            html.A(
                children='Surface Viewer',
                target= '_blank',
                href='https://brainbrowser.cbrain.mcgill.ca/surface-viewer#ct',
                style={'color': '#F012BE'}
            ),
            '.'
        ]
    )
], style={'margin': '0 auto'})

app.css.append_css({'external_url': 'https://codepen.io/plotly/pen/YeqjLb.css'})

@app.callback(
    Output('brain-graph', 'figure'),
    [Input('brain-graph', 'clickData'),
    Input('radio-options', 'value'),
    Input('colorscale-picker', 'colorscale')],
    [State('brain-graph', 'figure')])
def add_marker(clickData, val, colorscale, figure):

    if figure['data'][0]['name'] != val:
        if val == 'human':
            pts, tri=read_mniobj("realct.obj")
            intensities=np.loadtxt('realct.txt')
            figure['data']=plotly_triangular_mesh(pts, tri, intensities,
                                        colorscale=DEFAULT_COLORSCALE, flatshading=False,
                                        showscale=False, reversescale=False, plot_edges=False)
        elif val == 'human_atlas':
            pts, tri=read_mniobj("surf_reg_model_both.obj")
            intensities=np.loadtxt('aal_atlas.txt')
            figure['data']=plotly_triangular_mesh(pts, tri, intensities,
                                        colorscale=DEFAULT_COLORSCALE, flatshading=False,
                                        showscale=False, reversescale=False, plot_edges=False)
        elif val == 'mouse':
            pts, tri=read_mniobj("mouse_surf.obj")
            intensities=np.loadtxt('mouse_map.txt')
            figure['data']=plotly_triangular_mesh(pts, tri, intensities,
                                        colorscale=DEFAULT_COLORSCALE, flatshading=False,
                                        showscale=False, reversescale=False, plot_edges=False)
            pts, tri=read_mniobj("mouse_brain_outline.obj")
            outer_mesh = plotly_triangular_mesh(pts, tri)[0]
            outer_mesh['opacity'] = 0.5
            outer_mesh['colorscale'] = 'Greys'
            figure['data'].append(outer_mesh)
        figure['data'][0]['name'] = val

    elif clickData != None:
        if 'points' in clickData:
            marker = dict(
                x = [clickData['points'][0]['x']],
                y = [clickData['points'][0]['y']],
                z = [clickData['points'][0]['z']],
                mode = 'markers',
                marker = dict(size=15, line=dict(width=3)),
                name = 'Marker',
                type = 'scatter3d',
                text = ['Click point to remove annotation']
            )
            anno = dict(
                x = clickData['points'][0]['x'],
                y = clickData['points'][0]['y'],
                z = clickData['points'][0]['z'],
                font = dict(color = 'black'),
                bgcolor = 'white',
                borderpad = 5,
                bordercolor = 'black',
                borderwidth = 1,
                captureevents = True,
                ay = -50,
                arrowcolor = 'white',
                arrowwidth = 2,
                arrowhead = 0,
                text = 'Click here to annotate<br>(Click point to remove)',
            )
            if len(figure['data']) > 1:
                same_point_found = False
                for i, pt in enumerate(figure['data']):
                    if pt['x'] == marker['x'] and pt['y'] == marker['y'] and pt['z'] == marker['z']:
                        ANNO_TRACE_INDEX_OFFSET = 1
                        if val == 'mouse':
                            ANNO_TRACE_INDEX_OFFSET = 2
                        figure['data'].pop(i)
                        print('DEL. MARKER', i, figure['layout']['scene']['annotations'])
                        if len(figure['layout']['scene']['annotations']) >= (i-ANNO_TRACE_INDEX_OFFSET):
                            try:
                                figure['layout']['scene']['annotations'].pop(i-ANNO_TRACE_INDEX_OFFSET)
                            except:
                                pass
                        same_point_found = True
                        break
                if same_point_found == False:
                    figure['data'].append(marker)
                    figure['layout']['scene']['annotations'].append(anno)
            else:
                figure['data'].append(marker)
                figure['layout']['scene']['annotations'].append(anno)

    cs = []
    for i, rgb in enumerate(colorscale):
        cs.append([i/(len(colorscale)-1), rgb])
    figure['data'][0]['colorscale'] = cs

    return figure

@app.callback(
    Output('click-data', 'children'),
    [Input('brain-graph', 'clickData')])
def display_click_data(clickData):
    return json.dumps(clickData, indent=4)

@app.callback(
    Output('relayout-data', 'children'),
    [Input('brain-graph', 'relayoutData')])
def display_click_data(relayoutData):
    return json.dumps(relayoutData, indent=4)

if __name__ == '__main__':
    app.run_server(debug=True)
