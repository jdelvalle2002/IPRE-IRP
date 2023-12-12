import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import math

## implementaremos la matheuristicas ALNS

# utilizaremos networkx para representar el grafo, al que le aplicaremos un algoritmo
# que utilice ALNS para resolver un Inventory Routing Problem

# la idea es que el grafo represente la red de distribucion de una empresa, donde los nodos

demandas = random.sample(range(10, 100), 10)

SEED = 1234
random.seed(SEED)
# generamos el grafo de los clientes
G = nx.Graph()
G.add_nodes_from(range(10))
for i in range(10):
    for j in range(10):
        if i != j:
            decision = random.binomial(1, 0.7, seed=SEED)
            if decision == 1:
                G.add_edge(i, j, weight=random.randint(1, 10, seed=SEED))

# añadimos el nodo fuente

G.add_node(10)
# lo conectamos con la red mediante 3 aristas
utilizados = []
for l in range(3):
    k = random.randint(0, 9, seed=SEED)
    utilizados.append(k)
    G.add_edge(10, k, weight=random.randint(1, 10, seed=SEED))

# asignamos demandas a los nodos
for i in range(10):
    G.nodes[i]['demand'] = demandas[i]


