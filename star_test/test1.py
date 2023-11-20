import matplotlib.pyplot as plt
import networkx as nx

#from : https://stackoverflow.com/questions/64193989/how-do-i-curve-edges-in-networkx-graph

G = nx.Graph()

G.add_edge("Ted", "May", weight=0.5)
G.add_edge("Ted", "Ray", weight=5)
G.add_edge("Ted", "Chris", weight=1)
G.add_edge("Ted", "Sam", weight=20)
G.add_edge("Ted", "April", weight=8.8)
G.add_edge("Ted", "Ana", weight=0)


G.add_edge("Ana", "Ryan", weight=10)
G.add_edge("Ana", "Jim", weight=0.5)
G.add_edge("Ana", "Ben", weight=1)


pos = nx.circular_layout(G, scale=0.2)  # positions for all nodes
ax=plt.gca()

#for edge in G.edges():
#    source, target = edge
#    rad = -0.4
#    arrowprops=dict(lw=G.edges[(source,target)]['weight'],
#                    arrowstyle="-", 
#                    color='orange',
#                    connectionstyle=f"arc3,rad={rad}",
#                    linestyle= '-',
#                    alpha=0.3)
#    ax.annotate("",
#                xy=pos[source],
#                xytext=pos[target],
#                arrowprops=arrowprops
#               )
               
# https://stackoverflow.com/questions/25639169/networkx-change-color-width-according-to-edge-attributes-inconsistent-result
edges = G.edges()
#colors = [G[u][v]['color'] for u,v in edges]
weights = [G[u][v]['weight'] for u,v in edges]               
               
nx.draw_networkx_edges(
    G, pos, arrows=True, edge_color="red",  alpha=0.5 , width=weights,
    connectionstyle="arc3,rad=-0.2"  # <-- THIS IS IT
)
# nodes
nx.draw_networkx_nodes(G, pos, node_size=900, node_color="lightgray")
#sc.set_zorder(1) # .set_zorder(zorder)
nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")

plt.show()


############
## - with continuous color, but this will not be employed as would yield distinct color scales each time

# again
pos = nx.circular_layout(G, scale=0.2)  # positions for all nodes
ax=plt.gca()

colors_numbers = [2.2, 3, 4, 20, 5, 15, 6, 18,9]  # instead this, should retrieve the weiths into a list
nx.draw_networkx_nodes(G,pos,node_color='wheat', alpha=1,  node_size=700, edgecolors="black", node_shape="s")
colors_edges_l = ["lightgray", "royalblue", "pink", "salmon", "wheat", "blue", "darkcyan", "pink", "green"]
nx.draw_networkx_edges(G,pos, connectionstyle="arc3,rad=-0.3", alpha=1, arrows=True,
                               edge_color=colors_numbers,
                                        width=weights,                                  
                               edge_cmap=plt.cm.YlOrBr)  # use same **
# twice: this second one returns what can be read by colorbar
edges = nx.draw_networkx_edges(G,pos, 
                               edge_color=colors_numbers,
                                        width=0,                                  
                               edge_cmap=plt.cm.YlOrBr) # use same ** Purples ?
#edges                               
#nx.draw_networkx_edges(G,pos,  edge_color="salmon", alpha=0.3, arrows=True,
#                                        width=weights, connectionstyle="arc3,rad=-0.15")
plt.colorbar(edges)
plt.axis('off')
#plt.show()
plt.savefig("ugly.svg")
plt.close()
# problem: those with weak -log10padj are invisible! 
# because the color scale starts in white, and no scale starts with gray
# explore if I can do my own scale starting gray and ending brown


