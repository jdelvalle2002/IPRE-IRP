import numpy as np
import pandas as pd
import networkx as nx
import osmnx as ox
import math
import pandas as pd
import random
from matplotlib import pyplot as plt
import numpy as np
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from funciones import read_data, calcular_distancia, calcular_matriz_dist, calcular_largo_ruta 

random.seed(42)
np.random.seed(42)

ubis, cap_tpte, info_locales = read_data('IRP1.xlsx')

G = nx.DiGraph()
color_nodos = []
color_arcos = []
ancho_edges = []

for local in info_locales.itertuples():
    # print(local.i)
    G.add_node(f"N_{local.i}", Inv = local.I, Up = local.U, Low = local.L, Prod = local.r, h = local.h,
        coord_x = local.X, coord_y = local.Y, pos = (local.X, local.Y))
    if local.i != 0:
        color_nodos.append('blue') 
    else:
        color_nodos.append('red')

e=0
for local in G.nodes():
    for nodo in G.nodes():
        if local != nodo:
            decision = np.random.binomial(1, 0.7)
            if decision == 1 and local != 'N_0' and nodo != 'N_0':
                dist = calcular_distancia(G.nodes[local]['pos'], G.nodes[nodo]['pos'])
                G.add_edge(local, nodo, weight=dist)
                e +=1
                # ancho_edges.append(dist/1000)
                ancho_edges.append(0.25)
                color_arcos.append('gray')
            elif local == 'N_0':
                G.add_edge(local, nodo, weight=calcular_distancia(G.nodes[local]['pos'], G.nodes[nodo]['pos']))
                e +=1
                ancho_edges.append(1)
                color_arcos.append('black')
            else:
                ancho_edges.append(0)
    
plt.figure(figsize=(5,5))

pos=nx.get_node_attributes(G,'pos')
nx.draw(G, pos=pos, with_labels=True, node_size=18, font_size=15, node_color=color_nodos, width=ancho_edges, edge_color=color_arcos)
# plt.show()
print(G.nodes)
print(color_arcos, color_nodos, ancho_edges)