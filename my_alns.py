import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import math

## implementaremos la matheuristicas ALNS

# utilizaremos networkx para representar el grafo, al que le aplicaremos un algoritmo
# que utilice ALNS para resolver un Inventory Routing Problem

# la idea es que el grafo represente la red de distribucion de una empresa, donde los nodos

demandas = random.sample(range(1, 100), 10)

print(demandas)