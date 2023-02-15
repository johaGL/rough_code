import networkx as nx

def cyto2nx(nodes, edges):
    # ignore position
    G = nx.DiGraph()
    for dico in nodes:
        G.add_node(dico['data']['id'], score = dico['data']['score'])
    for dico in edges:
        x = dico['data']
        G.add_edge(x['source'], x['target'], edge_genes = x['edge_genes'])
    return G

def contractpath_nx(G, path):
    # get the genes of the path to be contracted
    genesset = set()
    # find the 'extremes' in the path, and the directionnality:
    gr = nx.DiGraph()
    gr.add_edges_from(path)
    startnodes = [n for n, d in gr.in_degree() if d == 0]
    endnodes = [n for n, d in gr.out_degree() if d == 0]
    #assert len(startnodes) == 1, "error, more than 1 startnodes"
    #assert len(endnodes) == 1, "error, more than 1 endnodes"
    if len(startnodes) == 1 and len(endnodes) == 1:
        try:
            for edge in path:
                tmpgdico = G.get_edge_data(edge[0], edge[1])   # {'edge_genes': [6, 7]}
                for k in tmpgdico['edge_genes']:
                    genesset.add(k)
            print(genesset)
            # connect end and start nodes in G and add all genes to this new shortcut
            F = G.copy()
            F.add_edge(startnodes[0], endnodes[0], edge_genes = list(genesset))
            # delete the path as demanded
            F.remove_edges_from(path)
            isol_ = list(nx.isolates(F))
            for i in isol_:
                #print('===>', i)
                F.remove_node(i)
            print("successfuly contracted path")
            return F

        except Exception as e:
            print("path could not be contracted, Error : \n", e)
            return G
    else:
        print("invalid paths")
        return G

   
def nx2cyto(F):
    nxnodes = nx.get_node_attributes(F, 'score')
    posx = 30
    posy = -100

    prepnodes = []
    for idnode in nxnodes:
        score = nxnodes[idnode]
        prepnodes.append((idnode, score, posx, posy))

    cytonodes = [ {
        'data': {'id': label, 'label': label},
        'position': {'x': 20 * lat, 'y': -20 * long}
    }
    for label, score, long, lat in prepnodes ]


    nxedges = F.edges()
    prepedges = []
    for ed in nxedges:
        tmpgdico = F.get_edge_data(ed[0], ed[1])
        genes = tmpgdico['edge_genes']
        prepedges.append((ed[0], ed[1], list(genes)))
    cytoedges = [
        {'data': {'source': source, 'target': target, 'edge_genes': edge_genes}}
        for source, target, edge_genes in prepedges
    ]

    return cytonodes, cytoedges




# note these templates to manipulate nodes and edges data
# as used by dash cytoscape :

#  node dictionaries
#  [{'data': {'id': 'Los Angeles', 'label': 'Los Angeles', 'score': '10'}, 'position': {'x': -2000, 'y': -600}},
#      {'data': {'id': 'Montreal', 'label': 'Montreal', 'score': '20'}, 'position': {'x': -2000, 'y': -600}},
#      {'data': {'id': 'Vancouver', 'label': 'Vancouver', 'score': '30'}, 'position': {'x': -2000, 'y': -600}},
#      {'data': {'id': 'Houston', 'label': 'Houston', 'score': '40'}, 'position': {'x': -2000, 'y': -600}},
#      {'data': {'id': 'Bordeaux', 'label': 'Bordeaux', 'score': '50'}, 'position': {'x': -2000, 'y': -600}}]
#
#
# edge dictionaries
#
#  [{'data': {'source': 'Los Angeles', 'target': 'Montreal', 'edge_genes': [1, 2]}},
#  {'data': {'source': 'Los Angeles', 'target': 'Vancouver', 'edge_genes': [3, 7]}},
#   {'data': {'source': 'Montreal', 'target': 'Houston', 'edge_genes': [9, 1]}},
#   {'data': {'source': 'Houston', 'target': 'Vancouver', 'edge_genes': [8, 5]}},
#    {'data': {'source': 'Vancouver', 'target': 'Houston', 'edge_genes': [11, 20]}},
#    {'data': {'source': 'Vancouver', 'target': 'Bordeaux', 'edge_genes': [4, 30]}},
    # {'data': {'source': 'Bordeaux', 'target': 'Los Angeles', 'edge_genes': [6, 7]}}]

if __name__ == '__main__':
    nodes = [{'data': {'id': 'Los Angeles', 'label': 'Los Angeles', 'score': '10'}, 'position': {'x': -2000, 'y': -600}},
     {'data': {'id': 'Montreal', 'label': 'Montreal', 'score': '20'}, 'position': {'x': -2000, 'y': -600}},
     {'data': {'id': 'Vancouver', 'label': 'Vancouver', 'score': '30'}, 'position': {'x': -2000, 'y': -600}},
     {'data': {'id': 'Houston', 'label': 'Houston', 'score': '40'}, 'position': {'x': -2000, 'y': -600}},
     {'data': {'id': 'Bordeaux', 'label': 'Bordeaux', 'score': '50'}, 'position': {'x': -2000, 'y': -600}}]


    edges = [{'data': {'source': 'Los Angeles', 'target': 'Montreal', 'edge_genes': [1, 2]}},
             {'data': {'source': 'Los Angeles', 'target': 'Vancouver', 'edge_genes': [3, 7]}},
             {'data': {'source': 'Montreal', 'target': 'Houston', 'edge_genes': [9, 1]}},
             {'data': {'source': 'Houston', 'target': 'Vancouver', 'edge_genes': [8, 5]}},
             {'data': {'source': 'Vancouver', 'target': 'Houston', 'edge_genes': [11, 20]}},
             {'data': {'source': 'Vancouver', 'target': 'Bordeaux', 'edge_genes': [4, 30]}},
             {'data': {'source': 'Bordeaux', 'target': 'Los Angeles', 'edge_genes': [6, 7]}}]

    G = cyto2nx(nodes,edges)
    print(nx2cyto(G))
    newnodes, newedges = nx2cyto(G)
    print(newnodes)
    print()
    print(newedges)

    mystery = [('Vancouver', 'Bordeaux'), ('Bordeaux', 'Los Angeles')]
    print(contractpath_nx(G, mystery))


# @app.callback(Output('avue', 'elements'),
#               Input('submit-val', 'n_clicks'),
#               State('cytoscape-selectedEdgeData-markdown', 'children'))
# def listen_edgesselection(n_clicks, children):
#     edgesstring = str(children)
#     if len(edgesstring) > 20:
#         strcli = edgesstring.split("Edges selection --> ")[1].strip()
#         edges_list = [tuple(k.split("%TO%")) for k in strcli.split(", ")]
#         G = cyto2nx(newnodes, newedges)
#         print(G)
#         print(edges_list)
#         F = contractpath_nx(G, edges_list)
#         print("====>" , F.nodes())
#         aka, oko = nx2cyto(F)
#         #newnodes = aka
#         #newedges = oko
#         return aka + oko
#     else:
#         print("nothing to contract")
#         return newnodes + newedges