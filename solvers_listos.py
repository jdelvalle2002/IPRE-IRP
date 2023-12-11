import numpy as np
import pandas as pd
import networkx as nx
import osmnx as ox
import math
import pandas as pd
import random
from matplotlib import pyplot as plt
import numpy as np

from funciones import (
    read_data,
    calcular_distancia,
    calcular_matriz_dist,
    calcular_largo_ruta,
    graficar_ruta,
)

# from instancia import ubis, cap_tpte, info_locales
# from instancia import G, color_nodos, color_arcos, ancho_edges
from instancia import *

random.seed(42)
np.random.seed(42)


def simular_demanda_diaria(data, dist="d", log=False):
    demanda = {}
    demanda[0] = 0
    for i in range(1, len(data)):
        if dist == "p":
            demanda[i] = int(round(np.random.poisson(data[i]), 0))
        elif dist == "n":
            demanda[i] = int(round(np.random.normal(data[i], data[i] / 10), 0))
        elif dist == "u":
            demanda[i] = int(round(np.random.uniform(data[i] / 2, data[i] * 1.5), 0))
        elif dist == "d":
            demanda[i] = int(round(data[i], 0))

    if log:
        print("Demanda del día: ", demanda)
        # print(i, ":", data[i])
    return demanda


def aplicar_politica_estacionaria(G, info_locales, matriz_dst):
    V = 0
    demanda = simular_demanda_diaria(info_locales["r"][1:], "d")
    for i in range(1, len(info_locales)):
        V += 0
    # Definimos el inventario del día de cada nodo


"""NEAREST NEIGHBOR ALGORITHM FOR VRP"""


def nearest_neighbor(
    G, dist_matrix, demands={}, capacity=float("inf"), disponibilidad=None
):
    """
    Apply the Nearest Neighbor heuristic to find initial routes for VRP.
    """
    # print('Ejecutando nearest_neighbor')
    if "N_0" not in demands.keys():
        demands = {f"N_{k}": v for k, v in demands.items()}
    nodos = []
    if disponibilidad is None:
        nodos = list(G.nodes)

    else:
        # recorremos la lista de nodos y nos quedamos con los que tienen disponibilidad
        for nodo in list(G.nodes):
            # print('¿Visito el nodo? ', nodo, disponibilidad[nodo])
            if disponibilidad[nodo]:
                nodos.append(nodo)
        # nodos = [nodo for nodo in list(G.nodes) if disponibilidad[nodo]]
    N = len(nodos)
    visitados = {nodo: False for nodo in nodos}
    rutas = []

    while sum(visitados.values()) < N:
        nodo_actual = "N_0"
        capacidad_actual = 0
        ruta = [nodo_actual]
        visitados[nodo_actual] = True
        # print('aca')
        while True:  # capacidad_actual + demands[nodo_actual] <= capacity:
            actual = ruta[-1]
            cercano = None
            min_dist = float("inf")
            # print('actual: ', actual)
            for vecino in G.neighbors(actual):
                # print('vecino: ', vecino)
                if vecino not in nodos or visitados[vecino]:
                    # print(f'el nodo {vecino} no está en la lista de nodos')
                    pass
                else:
                    # print(vecino, id_vecino, actual, dist_matrix[actual][id_vecino], demands[vecino])
                    # print(f"vecino: {id_vecino}, distancia: {dist_matrix[actual][id_vecino]} , min_dist: {min_dist}")
                    if (
                        dist_matrix[actual][vecino] < min_dist
                    ):  # demands[vecino] + capacidad_actual <= capacity
                        cercano = vecino
                        # print('cercano: ', cercano)
                        min_dist = dist_matrix[actual][vecino]
            # print(f'Nodo más cercano a {actual}: N_{cercano}')

            if cercano is None:
                # print('quibre con: ', ruta)
                break
            else:
                ruta.append(cercano)
                visitados[cercano] = True
                # capacidad_actual += demands[cercano] omitiremos la capacidad por ahora
        ruta.append("N_0")
        # print('Ruta actual: ', ruta, '\n')
        rutas.append(ruta)
    return rutas


"""TWO OPT ALGORITHM FOR VRP"""


def two_opt(ruta_inicial, matriz_dst, iters):
    ruta_2opt = ruta_inicial.copy()

    for k in range(iters):
        i, j = np.random.randint(1, len(ruta_2opt) - 1, size=2)
        if j < i:
            i, j = j, i

        nueva_ruta = ruta_2opt.copy()
        nueva_ruta[i:j] = ruta_2opt[j - 1 : i - 1 : -1]

        if calcular_largo_ruta(nueva_ruta, matriz_dst) < calcular_largo_ruta(
            ruta_2opt, matriz_dst
        ):
            # print("Ruta anterior: ", ruta_2opt)
            # print("Mejora encontrada -> Nueva ruta: ", nueva_ruta)
            ruta_2opt = nueva_ruta

    return ruta_2opt


def graficar_rutas(rutas, G):
    colores = ["blue", "green", "orange", "purple", "black", "red"]
    grafo = nx.DiGraph()
    nodos = list(G.nodes())
    color_nodos = []
    color_edges = []
    for nodo in nodos:
        grafo.add_node(nodo, pos=G.nodes()[nodo]["pos"])
        if nodo != "N_0":
            color_nodos.append("black")
        else:
            color_nodos.append("green")
    for ruta in rutas.values():
        if ruta != []:
            color = colores[-1]
            colores.remove(color)
            # print(f"Ruta: {ruta} en color: {color}")

            for i in range(len(ruta) - 1):
                n1, n2 = ruta[i], ruta[i + 1]
                grafo.add_edge(n1, n2, color=color)
                # print(f"Agregando arco: {n1} -> {n2} en color: {color}")
                color_edges.append(color)

    plt.figure(figsize=(5, 5))
    colors = nx.get_edge_attributes(grafo, "color").values()
    pos = nx.get_node_attributes(G, "pos")
    # edges = grafo.edges()
    nx.draw(
        grafo,
        pos=pos,
        with_labels=True,
        node_size=18,
        font_size=15,
        node_color=color_nodos,
        edge_color=colors,
    )
    plt.show()


def generar_ruta(G, matriz_dst, nodos_a_visitar):
    rutas_NN = nearest_neighbor(
        G, dist_matrix=matriz_dst, disponibilidad = nodos_a_visitar
    )
    # print( "Ruta NN: ")
    # print(rutas_NN[0])
    rutas_2opt = [two_opt(ruta, matriz_dst, 1000) for ruta in rutas_NN]
    # print( "Ruta 2O: ", rutas_2opt[0])
    # rutas = [rutas_NN[0], rutas_2opt[0]]
    # largo = calcular_largo_ruta(rutas[1], matriz_dst)
    rutas = rutas_2opt[0]
    largo = calcular_largo_ruta(rutas, matriz_dst)
    return rutas, largo  # retornamos la ruta 2opt y su largo

