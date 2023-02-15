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

hardedges = (
        ('la', 'mtl', [1, 2]),
        ('la', 'van', [3, 7]),
        ('mtl', 'hou', [9, 1]),
        ('hou', 'van', [8, 5]),
        ('van', 'hou', [11, 20]),
        ('van', 'bo', [4,30]),
        ('bo', 'la', [6,7])
    )

nodes = [
    {
        'data': {'id': short, 'label': label},
        'position': {'x': 20 * lat, 'y': -20 * long}
    }
    for short, label, long, lat in hardnodes
]

edges = [
    {'data': {'source': source, 'target': target, 'edge_genes' : edge_genes}}
    for source, target, edge_genes in hardedges
]

elements = nodes + edges

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
            'name': 'cose',
            'directed' : True
        }    ),

    html.P("a vue"),
    cyto.Cytoscape(
        id='avue',
        boxSelectionEnabled=True,  # shift + drag square to select enabled
        elements=newelems,
        style={'width': '100%', 'height': '200px'},
        stylesheet=default_stylesheet,
        layout={
            'name': 'cose',
            'directed' : True
        }),


    dcc.Markdown(id='cytoscape-selectedNodeData-markdown'),
    dcc.Markdown(id='cytoscape-selectedEdgeData-markdown'),
    html.H3(' ** '),
    # html.Pre(id='cytoscape-tapEdgeData-json', style=styles['pre']),
    # dcc.Markdown(id='cytoscape-selectedEdgeData-markdown'),

    html.Label('Text Input'),
    dcc.Input(id = 'true_false', value='F', type='text'),
    #html.Button('Submit', id='submit-val'),
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
              [Input('avue', 'selectedEdgeData')])
def displaylistofselnodes(data_list):
    if data_list is None:
        return "Edges selection -->"
    edges_list = [(data['source'], data['target']) for data in data_list]
    edges_liststr = ["%TO%".join(ktup) for ktup in edges_list ]
    #return "\n*".join(cities_list)
    return f'Edges selection -->:     {", ".join(edges_liststr)}'


# @app.callback(Output('cytoscape-tapEdgeData-json', 'children'),
#               [Input('avue', 'tapEdgeData')])
# def displayTapEdgeData(data):
#     return json.dumps(data, indent = 2)

# print(edges_list) # check if the variable is global ! == > no it si not !

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






if __name__ == '__main__':
    app.run_server(debug=True)