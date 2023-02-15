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



# @app.callback(Output('cytoscape-tapNodeData-json', 'children'),
#               [Input('cytoscape', 'tapNodeData')])
# def displayTapNodeData(data):
#     return json.dumps(data, indent=2)






# @app.callback(Output('avue', 'elements'),
#               [Input('cytoscape-selectedEdgeData-markdown', 'children'),
#                Input('submit-val', 'n_clicks')])
# def button_click(edgesString,  n_clicks):
#     if len(edgesString) >  20:
#         strcli = edgesString.split("Edges selection --> ")[1].strip()
#     else:
#         strcli = ''
#         print("=>"+strcli+"<==")
#     if strcli != '' and n_clicks == 1:
#         edges_list = [tuple(k.split("%TO%")) for k in strcli.split(", ")]
#         print(edges_list)
#         print(newedges)
#         print(newnodes)
#         return newnodes + newedges
#     else:
#         return newnodes + newedges

print(edgestocontract)

# @app.callback(Output('avue', 'elements'),
#               [Input('cytoscape-selectedEdgeData-markdown', 'children'),
#                Input('submit-val', 'n_clicks')])
# def listen_selectedEdgeData(edgesString):
#     if len(edgesString) > 20:
#         strcli = edgesString.split("Edges selection --> ")[1].strip()
#     else:
#         strcli = ''
#         print("=>" + strcli + "<==")
#         return newnodes + newedges
#
#     def button_click(strcli, n_clicks):
#         print(n_clicks)
#         edges_list = [tuple(k.split("%TO%")) for k in strcli.split(", ")]
#         print(edges_list)
#         print(newedges)
#         print(newnodes)
#         return newnodes + newedges





    #
    # if edges_list.split('-->')[1] == ' ' and n_clicks > 0 :
    #     print("doing nothing")
    #     n_clicks = 0
    # elif edges_list.split('-->')[1] == ' ' and n_clicks == 1:
    #     print(n_clicks)
    #     print(edges_list.split('-->')[1])
    #     return newnodes + newedges
    # else :
    #     print(n_clicks)
    #     print("more than one click is not allowed")

#
# @app.callback(Output('avue', 'elements'),
#                [Input('true_false', 'value')])
# def doweird(value):
#     print(value)
#     def givemecrazystuff():
#         newnodes = [
#             {
#                 'data': {'id': label, 'label': label},
#                 'position': {'x': 20 * lat, 'y': -20 * long}
#             }
#             for label, score, long, lat in (("X", 11, 30, -100), ("R", 20, 30, -100))
#         ]
#
#         newedges = [
#             {'data': {'source': source, 'target': target, 'edge_genes': edge_genes}}
#             for source, target, edge_genes in (('X', 'R', [1, 2]), ('R', 'X', [3, 4]))
#         ]
#         #print(newnodes + newedges)
#         return newnodes , newedges
#     if value == 'T':
#         newnodes , newedges = givemecrazystuff()
#         return newnodes + newedges
#     else:
#         return nodes + edges

