import numpy as np
import pandas as pd
import networkx as nx
import osmnx as ox
import math
import pandas as pd
import random
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm
from matplotlib.colors import ListedColormap
import funciones

def simular_demanda_previa(G, dist="n", T=100, ruido=0, d=1):
    """
    Función que simula la demanda previa de los locales.
    """
    # r = {nodo : nodo[1]['Prod'] for nodo in G.nodes(data=True)}
    if dist == "n":
        demandas = demanda_estable(G, T=T, ruido=ruido)
    elif dist == "c": # peak central
        demandas = demanda_peak_central(G, T=T, ruido=ruido)
    elif dist == "o": # oscilante
        demandas = demanda_oscilante(G, T=T, ruido=ruido, d=d)       
    elif dist == "d": # diagonal
        demandas = demanda_diagonal(G, T=T, ruido=ruido)    
    elif dist == "pir":
        demandas = demanda_piramide(G, T=T, ruido=ruido)
    return demandas

def demanda_estable(G, T=100, ruido=0):
    """
    Función que simula la demanda previa de los locales asumiendo que es constante.
    """
    g = G.copy()
    demandas = {nodo: [] for nodo in g.nodes() if nodo != "N_0"}
    # r = {nodo : nodo[1]['Prod'] for nodo in G.nodes(data=True)}
    for nodo in g.nodes(data=True):
        # print(nodo[0],nodo[1]['Prod'])
        if nodo[0] != "N_0":
            dem_pasadas = [
                max(
                    np.random.normal(
                        loc=nodo[1]["Prod"], scale=nodo[1]["Prod"] * 0.05
                    )
                    + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido),
                    0,
                )
                for _ in range(T)
            ]
            demandas[nodo[0]] = dem_pasadas
    return demandas
    
def demanda_peak_central(G, T=100, ruido=0):
    """
    Función que simula la demanda previa de los locales asumiendo que es constante, con un peak central por tramos.
    """
    g = G.copy()
    demandas = {nodo: [] for nodo in g.nodes() if nodo != "N_0"}
    # r = {nodo : nodo[1]['Prod'] for nodo in G.nodes(data=True)}
    for nodo in g.nodes(data=True):
        # print(nodo[0],nodo[1]['Prod'])
        if nodo[0] != "N_0":
            dem_pasadas = []
            for t in range(T):
                if 0.4*T <= t <= 0.45*T:
                    dem_pasadas.append(max(np.random.normal(loc=nodo[1]["Prod"]*1.25, scale=nodo[1]["Prod"] * 0.05) + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido), 0))
                elif 0.55*T <= t <= 0.6*T:
                    dem_pasadas.append(max(np.random.normal(loc=nodo[1]["Prod"]*1.25, scale=nodo[1]["Prod"] * 0.05) + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido), 0))
                elif 0.45*T <= t <= 0.55*T:
                    dem_pasadas.append(max(np.random.normal(loc=nodo[1]["Prod"]*1.5, scale=nodo[1]["Prod"] * 0.05) + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido), 0))
                else:
                    dem_pasadas.append(max(np.random.normal(loc=nodo[1]["Prod"], scale=nodo[1]["Prod"] * 0.05) + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido), 0))

            demandas[nodo[0]] = dem_pasadas
    return demandas

def demanda_oscilante(G, T=100, ruido=0, d=30):
    """
    Función que simula demanda que oscila en el tiempo con un comportamiento sinusoidal. P corresponde a la cantidad de peaks.
    """
    g = G.copy()
    demandas = {nodo: [] for nodo in g.nodes() if nodo != "N_0"}
    # r = {nodo : nodo[1]['Prod'] for nodo in G.nodes(data=True)}
    for nodo in g.nodes(data=True):
        if nodo[0] != "N_0":
            dem_pasadas = []
            for t in range(T):
                dem_pasadas.append(max(np.random.normal(loc=nodo[1]["Prod"], scale=nodo[1]["Prod"] * 0.05) * (1 + 0.5 * math.sin(math.pi * t * 2 / d)) + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido), 0))
            demandas[nodo[0]] = dem_pasadas
    return demandas

def demanda_diagonal(G, T=100, ruido = 0):
    g = G.copy()
    demandas = {nodo: [] for nodo in g.nodes() if nodo != "N_0"}
    # r = {nodo : nodo[1]['Prod'] for nodo in G.nodes(data=True)}
    for nodo in g.nodes(data=True):
        if nodo[0] != "N_0":
            dem_pasadas = []
            for t in range(T):
                dem_pasadas.append(max(np.random.normal(loc=nodo[1]["Prod"], scale=nodo[1]["Prod"] * 0.01) 
                                       * (0.1 * nodo[1]["Prod"] + 1.9*nodo[1]["Prod"] * t/T) # REVISAR IMPLEMENTACIÓN	
                                       + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido), 0))
            demandas[nodo[0]] = dem_pasadas
    return demandas

def demanda_piramide(G, T=100, ruido = 0):
    g = G.copy()
    demandas = {nodo: [] for nodo in g.nodes() if nodo != "N_0"}
    # r = {nodo : nodo[1]['Prod'] for nodo in G.nodes(data=True)}
    for nodo in g.nodes(data=True):
        if nodo[0] != "N_0":
            dem_pasadas = []
            for t in range(T):
                if t < T/2:
                    dem_pasadas.append(max(np.random.normal(loc=nodo[1]["Prod"], scale=nodo[1]["Prod"] * 0.01) 
                                       * (0.5 * nodo[1]["Prod"] + 1.5*nodo[1]["Prod"] * t/T) # REVISAR IMPLEMENTACIÓN	
                                       + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido), 0))
                else:
                    dem_pasadas.append(max(np.random.normal(loc=nodo[1]["Prod"], scale=nodo[1]["Prod"] * 0.01) 
                                       * (-0.5 * nodo[1]["Prod"] + 1.5*nodo[1]["Prod"] * t/T) # REVISAR IMPLEMENTACIÓN	
                                       + np.random.normal(loc=0, scale=nodo[1]["Prod"] * ruido), 0))
            demandas[nodo[0]] = dem_pasadas
    return demandas