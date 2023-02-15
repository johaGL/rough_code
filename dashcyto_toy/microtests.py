import dash
import dash_cytoscape as cyto
#import dash_html_components as html
# import dash_core_components as dcc
from dash import html
from dash import dcc
from pprint import pprint
from dash.dependencies import Input, Output, State
import json
from myfunctions import *

app = dash.Dash(__name__)

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

with open('data_nodes.txt', 'r') as f:
    h_n = f.readlines()
hardnodes = [k.strip().split('\t') for k in h_n]
posx = 30
posy = -100
hardnodes = [(i,j, posx, posy) for (i,j) in hardnodes]
#(hardnodes) # [('Los Angeles', 10, ...), ('Montreal', 20, ...), ...

hardedges = (
        ('Los Angeles', 'Montreal', [1, 2]),
        ('Los Angeles', 'Vancouver', [3, 7]),
        ('Montreal', 'Houston', [9, 1]),
        ('Houston', 'Vancouver', [8, 5]),
        ('Vancouver', 'Houston', [11, 20]),
        ('Vancouver', 'Bordeaux', [4,30]),
        ('Bordeaux', 'Los Angeles', [6,7])
    )

nodes = [
    {
        'data': {'id': label, 'label': label, 'score' : score},
        'position': {'x': 20 * lat, 'y': -20 * long}
    }
    for label, score, long, lat in hardnodes
]

edges = [
    {'data': {'source': source, 'target': target, 'edge_genes' : edge_genes}}
    for source, target, edge_genes in hardedges
]


default_stylesheet = [
    {
        'selector': 'node',
        'style': {
            #'background-color': '#FFFFFF',
            'label': 'data(label)',
            'color' : '#5A5A5A'
        }
    },

    {
        'selector' : 'edge',
        'style' : {
            'label' : 'data(edge_genes)',
            'curve-style' : 'bezier',
            'directed': True

        }
    }
]

 # initialize:

newedges = edges.copy()
newnodes = nodes.copy()
edgescontracted = []

print("global vars nodes edges")
print(len(newnodes))
print(len(newedges))

app.layout = html.Div([

    html.H2('Network handler'),
    html.P('the original network'),

    cyto.Cytoscape(
        id='cytoscape',
        boxSelectionEnabled=True,  # shift + drag square to select enabled
        elements=nodes + edges,
        style={'width': '100%', 'height': '250px'},
        stylesheet = default_stylesheet,
        layout={
            'name': 'cose',
            'directed' : True
        }    ),

    html.P(id='confirmationfield'),
    html.P("a vue"),
    cyto.Cytoscape(
        id='avue',
        boxSelectionEnabled=True,  # shift + drag square to select enabled
        elements=newnodes + newedges,
        style={'width': '100%', 'height': '250px'},
        stylesheet=default_stylesheet,
        layout={
            'name': 'cose',
            'directed' : True
        }),


    dcc.Markdown(id='cytoscape-selectedEdgeData-markdown'),
    html.Label('contract this path '),
    html.Button('Submit', id='submit-val'),
    dcc.Input(id='etc', value=None, type='text'),
    html.H4(' ** '),
    # html.Pre(id='cytoscape-tapEdgeData-json', style=styles['pre']),
    # dcc.Markdown(id='cytoscape-selectedEdgeData-markdown'),
    dcc.Markdown(id='cytoscape-selectedNodeData-markdown'),
    html.Label('Addthisnode'),
    dcc.Input(id = 'addthisnode', value='None', type='text'),

    # html.Label('Text Input'),
    # dcc.Input(id = 'true_false', value='F', type='text'),
    html.H4(' ** '),
    html.P('choose layout'),
    dcc.Dropdown(
            id = 'dropdown-update-layout',
            value = 'cose',
            clearable = False,
            options=[
                {'label': name.capitalize(), 'value': name}
                for name in ['grid', 'random', 'circle', 'cose', 'concentric']
            ]
        ),

    html.P('note : pending "bootstrap" css aesthetics'),
])

# https://dash.plotly.com/cytoscape/callbacks

@app.callback(Output('avue', 'layout'),
              [Input('dropdown-update-layout', 'value')])
def update_layout(layout):
    return {
        'name' : layout,
        'animate' : False
    }


@app.callback(Output('cytoscape-selectedNodeData-markdown', 'children'),
              [Input('avue', 'selectedNodeData')])
def displaylistofselnodes(data_list):
    if data_list is None:
        return "Nodes selection -->"
    cities_list = [data['label'] for data in data_list]
    #return "\n*".join(cities_list)
    return "Nodes selection -->:\n " + ", ".join(cities_list)


@app.callback(Output('cytoscape-selectedEdgeData-markdown', 'children'),
              [Input('avue', 'selectedEdgeData'), edgescontracted])
def displaylistofseledges(data_list, edgescontracted):
    if data_list is None:
        return "Edges selection -->"
    edges_list = [(data['source'], data['target']) for data in data_list]

    edges_liststr = ["%TO%".join(ktup) for ktup in edges_list ]
    #return "\n*".join(cities_list)
    edgesstring = f'Edges selection --> {", ".join(edges_liststr)}'
    return edgesstring



@app.callback(Output('avue', 'elements'),
              Input('submit-val', 'n_clicks'),
              State('cytoscape-selectedEdgeData-markdown', 'children'))
def listen_edgesselection(n_clicks, children):
    edgesstring = str(children)
    if len(edgesstring) > 20:
        strcli = edgesstring.split("Edges selection --> ")[1].strip()
        edges_list = [tuple(k.split("%TO%")) for k in strcli.split(", ")]
        G = cyto2nx(newnodes, newedges)
        print(G)
        print(edges_list)
        F = contractpath_nx(G, edges_list)
        print("====>" , F.nodes())
        aka, oko = nx2cyto(F)
        #newnodes = aka # not allowed
        #newedges = oko # not allowed
        return aka + oko
    else:
        print("nothing to contract")
        return newnodes + newedges



if __name__ == '__main__':
    app.run_server(debug=True)