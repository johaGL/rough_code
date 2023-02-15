import dash
import dash_cytoscape as cyto
#import dash_html_components as html
# import dash_core_components as dcc
from dash import html
from dash import dcc
from pprint import pprint
from dash.dependencies import Input, Output, State
import json

app = dash.Dash(__name__)

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

with open('data_nodes.txt', 'r') as f:
    h_n = f.readlines()
hardnodes = [i.strip().split('\t') for i in h_n]
hardnodes = [(i,j, 30, -100) for (i,j) in hardnodes]
print(hardnodes) # [('la', 'Los Angeles', 30, -100), ('mtl', 'Montreal', 30, -100), ('van', 'Vancouver', 30, -100), ('hou', 'Houston', 30, -100)]


nodes = [
    {
        'data': {'id': short, 'label': label},
        'position': {'x': 20 * lat, 'y': -20 * long}
    }
    for short, label, long, lat in hardnodes
]

edges = [
    {'data': {'source': source, 'target': target, 'edge_genes' : edge_genes}}
    for source, target, edge_genes in (
        ('mtl', 'la', [1, 2]),
        ('la', 'van', [3, 7]),
        ('mtl', 'van', [9, 1]),
        ('hou', 'van', [8, 5]),
        # ('mtl', 'la', 'faaa'),
        # ('la', 'van', 'k'),
        # ('mon', 'van', 'k'),
        # ('hou', 'van', 'k'),

    )
]

elements = nodes + edges

default_stylesheet = [
    {
        'selector': 'node',
        'style': {
            'background-color': '#BFD7B5',
            'label': 'data(label)'
        },

        'selector' : 'edge',
        'style' : {
            'label' : 'data(edge_genes)'
        }

    }
]

newelems = elements.copy() # initialize


app.layout = html.Div([

    html.H2('Network handler'),
    html.P('the original network'),

    cyto.Cytoscape(
        id='cytoscape',
        boxSelectionEnabled=True,  # shift + drag square to select enabled
        elements=elements,
        style={'width': '100%', 'height': '200px'},
        stylesheet = default_stylesheet,
        layout={
            'name': 'cose'
        }    ),

    html.P("a vue"),
    cyto.Cytoscape(
        id='avue',
        boxSelectionEnabled=True,  # shift + drag square to select enabled
        elements=newelems,
        style={'width': '100%', 'height': '200px'},
        stylesheet=default_stylesheet,
        layout={
            'name': 'cose'
        }),

    html.Pre(id='cytoscape-tapNodeData-json', style=styles['pre']),
    dcc.Markdown(id='cytoscape-selectedNodeData-markdown'),
    html.Pre(id='cytoscape-tapEdgeData-json', style=styles['pre']),
    dcc.Markdown(id='cytoscape-selectedEdgeData-markdown'),

    html.Label('Text Input'),
    dcc.Input(id = 'true_false', value='F', type='text'),
    html.Button('Submit', id='submit-val'),
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

    html.P('note : pending bootstrap aesthetics'),
])

# https://dash.plotly.com/cytoscape/callbacks

@app.callback(Output('avue', 'layout'),
              Input('dropdown-update-layout', 'value'))
def update_layout(layout):
    return {
        'name' : layout,
        'animate' : False
    }



@app.callback(Output('cytoscape-tapNodeData-json', 'children'),
              [Input('cytoscape', 'tapNodeData')])
def displayTapNodeData(data):
    return json.dumps(data, indent=2)

@app.callback(Output('cytoscape-tapEdgeData-json', 'children'),
              [Input('avue', 'tapEdgeData')])
def displayTapEdgeData(data):
    return json.dumps(data, indent = 2)

@app.callback(Output('avue', 'elements'),
               [Input('true_false', 'value')])
def doweird(value):
    print(value)
    def givemecrazystuff():
        newns = [
            {
                'data': {'id': short, 'label': label},
                'position': {'x': 20 * lat, 'y': -20 * long}
            }
            for short, label, long, lat in (("X", "ix", 30, -100), ("R", "ir", 30, -100))
        ]

        newedges = [
            {'data': {'source': source, 'target': target, 'edge_genes': edge_genes}}
            for source, target, edge_genes in (('X', 'R', [1, 2]), ('R', 'X', [3, 4]))
        ]
        print(newns + newedges)
        return newns + newedges
    if value == 'T':
        newelems = givemecrazystuff()
        return newelems
    else:
        return elements

#
# @app.callback(Output('avue', 'elements'),
#                [Input('submit-val', 'n_clicks')])
# def doweird(n_clicks):
#     def givemecrazystuff():
#         newns = [
#             {
#                 'data': {'id': short, 'label': label},
#                 'position': {'x': 20 * lat, 'y': -20 * long}
#             }
#             for short, label, long, lat in (("X", "ix", 30, -100), ("R", "ir", 30, -100))
#         ]
#
#         newedges = [
#             {'data': {'source': source, 'target': target, 'edge_genes': edge_genes}}
#             for source, target, edge_genes in (('X', 'R', [1, 2]), ('R', 'X', [3, 4]))
#         ]
#         print(newns + newedges)
#         return newns + newedges
#
#     newelems = givemecrazystuff()
#     print(n_clicks)
#     return newelems


if __name__ == '__main__':
    app.run_server(debug=True)