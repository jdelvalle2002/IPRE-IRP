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
from scipy.stats import norm

os.chdir(os.path.dirname(os.path.abspath(__file__)))


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

def crear_grafo_inicial(archivo="IRP1.xlsx", plot=False):
    ubis, cap_tpte, info_locales = read_data(archivo)

    G = nx.DiGraph()
    color_nodos = []
    color_arcos = []
    ancho_edges = []

    for local in info_locales.itertuples():
        G.add_node(
            f"N_{local.i}",
            Inv=local.I,
            Up=local.U,
            Low=local.L,
            Prod=local.r,
            h=local.h,
            coord_x=local.X,
            coord_y=local.Y,
            pos=(local.X, local.Y),
        )
        if local.i != 0:
            color_nodos.append("blue")
        else:
            color_nodos.append("red")

    for local in G.nodes():
        for nodo in G.nodes():
            if local != nodo:
                dist = calcular_distancia(G.nodes[local]["pos"], G.nodes[nodo]["pos"])
                G.add_edge(local, nodo, weight=dist)
                if local != "N_0" and nodo != "N_0":
                    ancho_edges.append(0.25)
                    color_arcos.append("gray")
                elif local == "N_0" or nodo == "N_0":
                    ancho_edges.append(1)
                    color_arcos.append("black")
    if plot:
        plt.figure(figsize=(5, 5))
        pos = nx.get_node_attributes(G, "pos")
        nx.draw(
            G,
            pos=pos,
            with_labels=True,
            node_size=18,
            font_size=15,
            node_color=color_nodos,
            width=ancho_edges,
            edge_color=color_arcos,
        )
        plt.show()

    return G, ubis, cap_tpte, info_locales

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
        sum_dist += matriz_dist[nodo][next_nodo]

    return sum_dist

def graficar_ruta(ruta, G):
    grafo = nx.DiGraph()
    nodos = list(G.nodes())

    color_nodos = []
    for nodo in nodos:
        
        if nodo in ruta and nodo != 'N_0':
            grafo.add_node(nodo, pos=G.nodes()[nodo]['pos'])
            color_nodos.append('red')
        elif nodo == 'N_0':
            grafo.add_node(nodo, pos=G.nodes()[nodo]['pos'])
            color_nodos.append('green')
        else:
            grafo.add_node(nodo, pos=G.nodes()[nodo]['pos'])
            color_nodos.append('gray')

    for i in range(len(ruta)-1):
        n1, n2 = ruta[i], ruta[i+1]

        grafo.add_edge(n1, n2, color='red')

    plt.figure(figsize=(5,5))
    pos=nx.get_node_attributes(grafo,'pos')
    nx.draw(grafo, pos=pos, with_labels=True, node_size=18, font_size=15,  node_color = color_nodos)
    plt.show()

def graficar_rutas(rutas, G):
    grafo = nx.DiGraph()
    nodos = list(G.nodes())
    color_nodos = []
    for nodo in nodos:
        id_nodo = int(nodo[2:])
        grafo.add_node(nodo, pos=G.nodes()[nodo]['pos'])
        if id_nodo != 0:
            color_nodos.append('red')
        else:
            color_nodos.append('green')
    for ruta in rutas:
        for i in range(len(ruta)-1):
            n1, n2 = ruta[i], ruta[i+1]
            grafo.add_edge(f'N_{n1}', f'N_{n2}', color='red')

    plt.figure(figsize=(5,5))
    pos=nx.get_node_attributes(G,'pos')
    nx.draw(grafo, pos=pos, with_labels=True, node_size=18, font_size=15,  node_color = color_nodos)
    plt.show()

def calcular_matriz_dist_alns(G):
    matriz_distancias = {int(nodo[2:]): {} for nodo in G.nodes()}
    for arco in G.edges(data=True):
        n1, n2, data = arco
        n1 = int(n1.split("_")[-1])
        n2 = int(n2.split("_")[-1])
        matriz_distancias[n1][n2] = data['weight']

    for n1 in G.nodes():
        n1 = int(n1.split("_")[-1])
        for n2 in G.nodes():
            n2 = int(n2.split("_")[-1])
            if n2 not in matriz_distancias[n1].keys():
                matriz_distancias[n1][n2] = np.inf
    return matriz_distancias

def simular_demanda_previa(G, dist = 'n', T =100, ruido = 0):
    """
    Función que simula la demanda previa de los locales.
    """
    g = G.copy()
    demandas = {nodo : [] for nodo in g.nodes() if nodo != 'N_0'}
    # r = {nodo : nodo[1]['Prod'] for nodo in G.nodes(data=True)}
    if dist == 'n':
        for nodo in g.nodes(data=True):
            # print(nodo[0],nodo[1]['Prod'])
            if nodo[0] != 'N_0':
                dem_pasadas = [max(
                    np.random.normal(loc = nodo[1]['Prod'], scale = nodo[1]['Prod'] * 0.05) 
                    + np.random.normal(loc = 0, scale = nodo[1]['Prod'] * ruido)
                    ,0) for _ in range(T)]
                demandas[nodo[0]] = dem_pasadas
    return demandas

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def SEDA(datos, historia = False, alpha=0.1, beta=0.1, theta=0.5):
    '''
    Suavizamiento Exponencial Doble Amortiguado
    Aplica el método de suavizamiento exponencial doble a una serie de datos,
    específicmamente el Método de Holt Damped.
    '''
    I = [datos[0]]
    S = [datos[1]-datos[0]]
    for i in range(1,len(datos)):
        I.append(alpha*datos[i] + (1-alpha)*(I[i-1]+ theta*S[i-1]))
        S.append(beta*(I[i]-I[i-1])+(1-beta)*S[i-1])
    
    y = I[-1] + theta * S[-1]
    if historia == False:
        return y
    elif historia == True:
        I.append(y)
        return I
    elif historia == 'S':
        return I, S
    
def pronostico_SEDA(datos, T, pron = False, alpha=0.1, beta=0.1, theta=0.5):
    '''
    Devuelve un pronóstico para los siguientes T periodos mediante Suavizamiento Exponencial Doble Amortiguado
    '''
    I,S = SEDA(datos, historia = 'S', alpha=alpha, beta=beta, theta=theta)
    pronostico = []
    for i in range(T):
        y = I[-1] + theta * S[-1]
        pronostico.append(y)
        I.append(alpha*y + (1-alpha)*(I[-1]+ theta*S[-1]))
        S.append(beta*(I[-1]-I[-2]) + (1-beta)*S[-1])
    
    return pronostico

def IC_nrm(mu, sd, M=1000, alfa = 0.95):
    """
    Función que calcula el intervalo de confianza para una distribución normal.
    """
    limite_inferior = mu + norm.ppf((1 - alfa) / 2) * sd / math.sqrt(M)
    limite_superior = mu - norm.ppf((1 - alfa) / 2) * sd / math.sqrt(M)
    
    return limite_inferior, limite_superior

def ejecutar_ruta(G,ruta,matriz_dst):
    """
    Función que simula la ejecución de una ruta.
    """
    g = G.copy()
    ruta = ruta.copy()
    ruta.pop(0)
    ruta.pop(-1)
    ruta = [int(nodo[2:]) for nodo in ruta]
    # distancia = calcular_largo_ruta(ruta, matriz_dst)
    stock = 0
    for nodo in ruta:
        stock += g.nodes[f'N_{nodo}']['Up'] - G.nodes[f'N_{nodo}']['Inv']
        g.nodes[f'N_{nodo}']['Inv'] = G.nodes[f'N_{nodo}']['Up']
    g.nodes['N_0']['Inv'] -= stock
    return g, stock

def realizacion_demanda(G0, ruido = 0.05):
    """
    Función que simula la demanda de los locales para un determinado periodo.
    """
    grafo = G0.copy()
    demandas = {nodo : [] for nodo in grafo.nodes() if nodo != 'N_0'}
    insatisfecho = 0
    for nodo in grafo.nodes(data=True):
        # print(nodo)
        if nodo[0] != 'N_0':
            dem = max(
                np.random.normal(loc = nodo[1]['Prod'], scale = nodo[1]['Prod'] * 0.05) 
                + np.random.normal(loc = 0, scale = nodo[1]['Prod'] * 0.05)
                ,0)
            demandas[nodo[0]] = dem
            if dem <= grafo.nodes[nodo[0]]['Inv']:
                grafo.nodes[nodo[0]]['Inv'] -= dem

            else:
                grafo.nodes[nodo[0]]['Inv'] = 0
                insatisfecho += dem - grafo.nodes[nodo[0]]['Inv']
        # print(nodo[0],nodo[1]['Inv'])
    # print(demandas)
    return grafo, demandas, insatisfecho

def reaccion_inventario(graf, mu, sd, alfa = 0.05):
    """
    Función que verifica que locales deben ser visitados en base a su inventario actual. 
    En caso de que el inventario se encuentre bajo el umbral de tolerancia, se retorna True.
    """
    grafo = graf.copy()
    visitas = {nodo : False for nodo in G.nodes()}
    for nodo in grafo.nodes(data=True):
        id_nodo = int(nodo[0][2:])-1
        media = mu[id_nodo]
        desviacion = sd[id_nodo]
        s = media + norm.ppf((1 - alfa)/2)* desviacion  #Stock de seguridad
        # print(f'{nodo[1]["Inv"]}, s{int(nodo[0][2:])} = {s}, {nodo[1]["Inv"] <= s}')
        if nodo[1]["Inv"] <= s:
            visitas[nodo[0]] = True
            # print(f'Visitar {nodo[0]}')
    # print(visitas)
    return visitas


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