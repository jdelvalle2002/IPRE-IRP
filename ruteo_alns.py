from alns import ALNS
from alns.accept import RecordToRecordTravel, SimulatedAnnealing
from alns.select import RouletteWheel
from alns.stop import MaxRuntime
import numpy as np
import pandas as pd
import networkx as nx
import osmnx as ox
import math
import pandas as pd
import random
from matplotlib import pyplot as plt
import numpy as np

from instancia import crear_grafo_inicial
from solvers_listos import (simular_demanda_diaria)

from funciones import (
    read_data,
    calcular_distancia,
    calcular_matriz_dist,
    calcular_largo_ruta,
    graficar_ruta,
    calcular_matriz_dist_alns
)

random.seed(42)
np.random.seed(42)

# Cargar datos
G, ubis, cap_tpte, info_locales = crear_grafo_inicial(archivo= 'IRP1.xlsx' ,plot=False)
matriz_dst = calcular_matriz_dist_alns(G)
demandas = simular_demanda_diaria(list(info_locales["r"]), dist="n")
print(demandas)

degree_of_destruction = 0.2
nodes_to_destroy = int(np.ceil(degree_of_destruction * len(G.nodes)))


'''
Definición de la clase IRPState
'''
class IRPState:
    """
    Estado de solución del IRP
    """

    def __init__(self, routes, unnasigned=None, matriz_dst = matriz_dst):
        self.routes = routes  # será solo una ruta ya que es solo un vehículo
        self.unnasigned = unnasigned if unnasigned is not None else []

    def copy(self):
        return IRPState(self.routes.copy(), self.unnasigned.copy())

    def objective(self):
        ruta = self.routes[0][:]
        ruta = [int(x.split("_")[-1]) for x in ruta]
        return calcular_largo_ruta(ruta, matriz_dst)

    @property
    def cost(self):
        return self.objective()

    def find_route(self, node):
        for route in self.routes:
            if node in route:
                return route

'''
Operaciones y funciones auxiliares
'''
def random_removal(state, random_state):
    destroyed = state.copy()

    for node in random_state.choice(
        range(1, len(G.nodes)), nodes_to_destroy, replace=False
    ):
        node = f"N_{node}"
        if node == "N_0":
            print("No se puede eliminar el depósito")
        destroyed.unnasigned.append(node)
        route = destroyed.find_route(node)
        if route is not None:
            route.remove(node)
    
    non_empty_routes = [ruta for ruta in destroyed.routes if len(ruta) > 0]

    destroyed.routes = non_empty_routes
    return destroyed

def greedy_repair(state, random_state):
    repaired = state.copy()

    random_state.shuffle(repaired.unnasigned)

    while len(repaired.unnasigned) > 0:
        node = repaired.unnasigned.pop()
        route, idx = best_insert(node, repaired)
        if route is not None:
            route.insert(idx, node)
        else:
            repaired.routes.append([node])
    return repaired

def best_insert(node, state):
    best_cost, best_route, best_idx = float("inf"), None, None
    for route in state.routes:
        for i in range(1, len(route) + 1):
            if can_insert(node, route):
                cost = insertion_cost(node, route, i)

                if cost < best_cost:
                    best_cost = cost
                    best_route = route
                    best_idx = i

    return best_route, best_idx

def can_insert(node, route, demanda, cap = 871):
    total = sum([demanda[int(nodo[2:])] for nodo in route])
    return total + demanda[int(node[2:])] <= cap   ####OJO ACÁ CON EL CAP, SI QUIERO CAMBIAR LA CAPACIDAD DEL VEHÍCULO TENGO QUE MODIFICARLO ACÁ

def insertion_cost(node, route, idx):
    pred = "N_0" if idx == 0 else route[idx - 1]
    succ = "N_0" if idx == len(route) else route[idx]
    pred = int(pred.split("_")[-1])
    succ = int(succ.split("_")[-1])
    node = int(node.split("_")[-1])
    dist = matriz_dst[pred][node] + matriz_dst[node][succ] - matriz_dst[pred][succ]
    return dist

'''
Simulación de demanda para un día
'''
def simular_demanda_diaria(data, dist="d", log=False):
    demanda = {}
    demanda['N_0'] = 0
    for i in range(1, len(data)):
        if dist == "p":
            demanda[f'N_{i}'] = int(round(np.random.poisson(data[i]), 0))
        elif dist == "n":
            demanda[f'N_{i}'] = int(round(np.random.normal(data[i], data[i] / 10), 0))
        elif dist == "u":
            demanda[f'N_{i}'] = int(round(np.random.uniform(data[i] / 2, data[i] * 1.5), 0))
        elif dist == "d":
            demanda[f'N_{i}'] = int(round(data[i], 0))

    if log:
        print("Demanda del día: ", demanda)
        # print(i, ":", data[i])
    return demanda

def neighborhood(G, node, distancias):
    # ordenaremos los nodo por distancia al nodo actual
    return sorted(G.neighbors(node), key=lambda x: distancias[int(node[2:])][int(x[2:])])


def nearest_neighbor_adapted(grafo, demandas, distancias, cap):
    """
    Adaptación del algoritmo Nearest Neighbor descrito anteriormente para que se utilice con la
    clase IRPState
    """
    # print(f'Capacidad del vehículo: {cap} | demanda total = {sum(demandas.values())}')
    G0 = grafo.copy()
    routes = []
    available = list(G0.nodes)
    available.remove("N_0")
    route = ["N_0"]  # siempre empezamos en el nodo 0
    route_demands = 0
    distancias = calcular_matriz_dist_alns(G0)

    while available:
        current = route[-1]

        my_neighbors = neighborhood(G0, current, distancias=distancias)
        while my_neighbors[0] in route:
            my_neighbors = my_neighbors[1:]
        nearest = my_neighbors[0]
        
        # print(f'uso de vehículo: {route_demands} | demanda del nodo {nearest}: {demandas[int(nearest[2:])]}')
        if route_demands + demandas[int(nearest[2:])] <= cap:
            route.append(nearest)
            route_demands += demandas[int(nearest[2:])]
            # print('nodesd:', nodes, 'nearest:', nearest)
            available.remove(nearest)
        else:
            break
    routes.append(route)
    return IRPState(routes)

def ruteo_ALNS(grafo, demandas, cap, F=1, ruta = False, MRT = 1):
    matriz_dst = calcular_matriz_dist_alns(grafo)
    init = nearest_neighbor_adapted(grafo, demandas, distancias = matriz_dst, cap=cap)

    alns = ALNS(np.random.RandomState(42))
    alns.add_destroy_operator(random_removal)
    alns.add_repair_operator(greedy_repair)

    select = RouletteWheel([10, 5, 2, 0], 0.7, 1, 1)
    # accept = RecordToRecordTravel.autofit(init.objective(), 0.02, 0, 9000)
    accept = SimulatedAnnealing(
        start_temperature=20000,
        end_temperature=0.1,
        step=0.9993,
        method="exponential",
    )
    stop = MaxRuntime(MRT)

    result = alns.iterate(init, select, accept, stop)
    init = result.best_state
    if ruta:
        solution = result.best_state
        route = solution.routes[0]
        return route
    
    else:
        return result

if __name__ == "__main__":
    result = ruteo_ALNS(G, demandas, cap=np.inf, F=1)
    solution = result.best_state
    objective = solution.objective()
    # pct_diff = 100 * (objective - bks.cost) / bks.cost

    print(f"Best heuristic objective is {objective}.")
    for idx, route in enumerate(solution.routes):
        print(f"Route {idx}:", route + ['N_0'])
        graficar_ruta(G= G, ruta= route + ['N_0'])