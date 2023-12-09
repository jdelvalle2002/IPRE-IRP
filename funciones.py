## Codigo basado en el trabajo de https://medium.com/@writingforara/solving-vehicle-routing-problems-with-python-heuristics-algorithm-2cc57fe7079c
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

random.seed(42)
np.random.seed(42)


def read_data(filename):
    """
    Read coordinates and demand values from a specific sheet in an Excel file.
    Assumes the data is in columns labeled 'X', 'Y', and 'Demand'.
    """
    df_grl = pd.read_excel(filename, sheet_name='GRL') #N	H	CP
    df_locales = pd.read_excel(filename, sheet_name='RTL') #i	X	Y	I	U	L	r	h

    ubicaciones = df_locales[['X', 'Y']].values 
    ubis = [list(ubicaciones[i]) for i in range(len(ubicaciones))]
    # info_locales = df_locales[['i','I','U','L','r','h']].values
    cap_tpte = df_grl['CP'].values[0]

    # return df_locales
    return ubis, cap_tpte, df_locales

def calcular_distancia(n1, n2):

    x1, y1 = n1[0], n1[1]
    x2, y2 = n2[0], n2[1]

    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def calcular_matriz_dist(G):
    matriz_distancias = {nodo: {} for nodo in G.nodes()}
    for arco in G.edges(data=True):
        n1, n2, data = arco
        matriz_distancias[n1][n2] = data['weight']

    for n1 in G.nodes():
        for n2 in G.nodes():
            if n2 not in matriz_distancias[n1].keys():
                matriz_distancias[n1][n2] = np.inf
    return matriz_distancias

def calcular_largo_ruta(ruta, matriz_dist):
    """
    Calculate the total distance of a given ruta using the distance matrix.
    """
    sum_dist = 0
    num_puntos = len(ruta)

    for i in range(num_puntos - 1):
        nodo = ruta[i]
        next_nodo = ruta[i + 1]
        sum_dist += matriz_dist[nodo, next_nodo]

    return sum_dist

# ubis, cap_tpte, info_locales = read_data('IRP1.xlsx')

# G = nx.DiGraph()
# color_nodos = []
# color_arcos = []
# ancho_edges = []

# for local in info_locales.itertuples():
#     # print(local.i)
#     G.add_node(f"N_{local.i}", Inv = local.I, Up = local.U, Low = local.L, Prod = local.r, h = local.h,
#         coord_x = local.X, coord_y = local.Y, pos = (local.X, local.Y))
#     if local.i != 0:
#         color_nodos.append('blue') 
#     else:
#         color_nodos.append('red')

# e=0
# for local in G.nodes():
#     for nodo in G.nodes():
#         if local != nodo:
#             decision = np.random.binomial(1, 0.7)
#             if decision == 1 and local != 'N_0' and nodo != 'N_0':
#                 dist = calcular_distancia(G.nodes[local]['pos'], G.nodes[nodo]['pos'])
#                 G.add_edge(local, nodo, weight=dist)
#                 e +=1
#                 # ancho_edges.append(dist/1000)
#                 ancho_edges.append(0.25)
#                 color_arcos.append('gray')
#             elif local == 'N_0':
#                 G.add_edge(local, nodo, weight=calcular_distancia(G.nodes[local]['pos'], G.nodes[nodo]['pos']))
#                 e +=1
#                 ancho_edges.append(1)
#                 color_arcos.append('black')
#             else:
#                 ancho_edges.append(0)
    
# #     for j in range(len(ubis)):
# #         if i != j:
# #             decision = np.random.binomial(1, 0.7)
# #             if decision == 1:
# #                 print(i, j)
# #                 G.add_edge(i, j, weight=random.randint(1, 10))

# plt.figure(figsize=(5,5))
# # nx.draw(G, with_labels=True, node_size=18, node_color=color_nodos, font_size=15, edge_color=ancho_edges, width=ancho_edges, edge_cmap=plt.cm.Greys)
# pos=nx.get_node_attributes(G,'pos')
# nx.draw(G, pos=pos, with_labels=True, node_size=18, font_size=15, node_color=color_nodos, width=ancho_edges, edge_color=color_arcos)
# plt.show()