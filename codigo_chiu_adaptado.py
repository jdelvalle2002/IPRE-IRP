# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 15:27:41 2022

@author: 56945
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 20:28:10 2021

@author: 56945

"""
#%%

"""
Les envío el código comentado. Hay dos scripts:

    Exp 6 Generador de instancias: se generan instancias con demanda Poisson
    Cluster6_poisson: lee cada instanciia generada anteriormente, corre la simulacion y guarda los resultados en un archivo .csv.

Para correr uns simulación con un modelo nuevo, hay que modificar el scrip Cluster6_poisson en: # Modificacion 1, # Modificacion 2 y # Modificacion 3.

    # Modificacion 1: Se definen en la simulación los modelos que se pueden correr
    # Modificacion 2: Se define el modelo a resolver y su método de solución
    # Modificacion 3: Se señala a la simulación que se corra X método de solución

Consideraciones generales:

    Para que corra se necesita Gurobi
    Está la simulación con una demanda Poisson
    Las demandas de la simulacion son generadas en Exp 6 Generador de instancias.py
    Todos los sets corresponden a diccionarios. Para ver su estructura conviene revisar Exp 6 Generador de instancias.py
    La simulación considera que estamos en un período t_0 y planifica para los períodos t_0 + 1, ..., t_0 + j. Donde en todos los períodos para los que se planifica la demanda es desconocida (se planifica utilizando la media de la distribución de probabilidad de la demanda).
    El código de la simulación consta de 3 grandes partes:
        Definición de funciones para la toma de decision
        Definicion de objetos claves:
            Objeto IRP: modela el estado del sistema y las transiciones
            Objeto Simulación: solamente corre la simulación
        Leer las instancias y correr la simulación 
"""

##################################################################################
########################### PASO 1: IMPORTAR LIBRERÍAS ###########################
##################################################################################

import imageio
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import random
import statistics
import time
import sys
import math

from gurobipy import *
from statistics import mean
<<<<<<< HEAD
=======
from funciones import *
>>>>>>> 7958b871f0e863ac9f6a46fd3d57cc9f2c2e5468




def inverse_function(lambda_):
    sample =  np.random.poisson(lambda_, 100000)
    value = np.percentile(sample, 95)
    return value


#######################################################################################################
########################### PASO 5: DEFINIR DEFINIR OBJETO DE SIMULACION P1 ###########################
#######################################################################################################

class IRP:
    def __init__(self, pickle_name, real_demands, option = None, cost_params = None):

        self.t = 0
        unpickled = pd.read_pickle(pickle_name)

        #Sets
        self.T = unpickled['Sets']['T'] #periodos
        self.C = unpickled['Sets']['C'] #clientes
        self.N = unpickled['Sets']['N'] #nodos
        self.A = unpickled['Sets']['A'] #arcos
        self.K = unpickled['Sets']['K'] #vehiculos

        #Parameters
        h = unpickled['Param']['h'] #costo de inventario
        self.Q = unpickled['Param']['Q'] #capacidad de vehiculo
        self.cap = unpickled['Param']['cap'] #capacidad de bodega
        self.pos = unpickled['Pos']['pos'] #matriz de posiciones en el plano cartesiano
        self.c = unpickled['Pos']['cost'] #matriz de costos de transporte
        self.h = {i: h[i] for i in self.N for t in self.T} #costo de inventario
        self.a = unpickled['Param']['a'] #costo de demanda insatisfecha
        self.r = unpickled['Param']['r'] #cantidad de unidades que llegan al depot en cada periodo

        #Demands
        self.mean = unpickled['Probs']['mu'] #demanda promedio (asume que saca la funcion inversa para obtener los stock de seguridad)
        #self.deviation = unpickled['Probs']['sigma']

        #safety_stock
        self.lambda_dict_original = unpickled['Probs']['mu']
        self.inverse_dict = {}
        self.lambda_dict = {}

        for customer in self.lambda_dict_original:
            for t in [1, 2, 3]:
                self.lambda_dict[(customer, t)] = t*self.lambda_dict_original[customer] #se asume demanda poisson

        self.safety_stock = {(i,c): inverse_function(self.lambda_dict[i, c]) - c*self.mean[i] for i in self.lambda_dict_original for c in (1,2,3)} #determinar stock de seguridad

        #Status
        self.status = {i: unpickled['Param']['init'][i] for i in self.N} #estado del sistema

        #real_data
        self.real_demands = real_demands #demandas simuladas a priori

        self.real_demands = {t: {i: int(self.real_demands[t][i][0]) for i in self.real_demands[t].keys()} for t in self.real_demands.keys()}

#        print("self.real_demands = real_demands")
#        print(f"real_demands= {self.real_demands}")

        #Revealed demands
        self.reveal = {i: 0 for i in self.C} #demanda revelada en cada periodo

        #Statistics
        self.inventory = {}         #se guardan los inventarios
        self.lost_sales = {}        #se guardan la dda insatisfecha
        self.routing_cost = {}      #se guardan los costos de ruteo
        self.final_stats = {}

        #Actions
        self.x_ = {}                #se guarda decision del modelo
        self.q_ = {}                #se guarda decision del modelo
        self.y_ = {}                #se guarda decision del modelo

        #History #decisiones historicas
        self.inventory_hist = {}
        self.lost_sales_hist = {}
        self.routes_hist = {}
        self.quantity_hist = {}
        self.sigma_hist = {} #flujo
        self.visits_hist = {}
        self.pred_hist = {}
        self.fill_rate_hist = {}
        self.gap_hist = {}

        #Suma de costos
        self.inventory_costs = 0
        self.lost_sales_costs = 0
        self.routing_costs = 0

        #Option
        self.option = option

        #Cost function
        self.cost_function = cost_params


    ###############################################################
    ############# Acá se toca para correr la simulacion ###########
    ###############################################################

    # Modificacion 1
    def action(self):

        sets = {'C': self.C, 'N': self.N, 'A': self.A, 'K': self.K}
        params = {'Q': self.Q, 'cap': self.cap, 'c': self.c, 'h': self.h, 'a': self.a, 'r': self.r, 'ss': self.safety_stock}
        probs = {'mu': self.mean}

        status = self.status

        ###REVISAR OPCIONES

        if self.option == 5:
            self.x_, self.q_, self.y_, gap, pred = Solve_IRP_SS(sets, params, probs, status, t0 = self.t, horizon = 3)


        self.gap_hist[self.t + 1] = gap
        self.visits_hist[self.t + 1] = self.y_
        self.pred_hist[self.t + 2] = pred
        pass

    # se revela la demanda
    def reveal_demands(self):

#        print("\n--------------REVEALING DEMANDS")

        self.reveal = self.real_demands[self.t + 1]

#        print("self.reveal = self.real_demands[self.t + 1]")
#        print(f"self.real_demands= {self.real_demands}")
#        print(f"self.reveal= {self.reveal}")
#        print(f"self.real_demands= {self.real_demands}")

        pass

    # actualizacion del estado
    def update_status(self):

#        print("\n--------------ACTUALIZANDO STATUS")

        status = {i: np.round(max(self.status[i] + sum(self.q_[i, k] for k in self.K) - self.reveal[i], 0), 4) for i in self.C}
        auxstatus = {i: np.round(max(self.status[i] + sum(self.q_[i, k] for k in self.K) - self.reveal[i], 0), 4) for i in self.C}

#        print("status = {i: np.round(max(self.status[i] + sum(self.q_[i, k] for k in self.K) - self.reveal[i], 0), 4) for i in self.C}")

#        print(f"self.q_= {self.q_}")
#        print(f"self.status= {self.status}")
#        print(f"self.K= {self.K}")
#        print(f"self.reveal= {self.reveal}")
#        print(f"self.C= {self.C}")

#        print(f"status inicial= {status}")
#        print(f"aux status inicial (sin np.round)= {auxstatus}")

        status[0] = np.round(max(self.status[0] + self.r - sum(self.q_[i, k] for k in self.K for i in self.C), 0), 4)

#        print(f"status 0 = {status}")

        I = {i: np.round(max(self.status[i] + sum(self.q_[i, k] for k in self.K) - self.reveal[i], 0), 4) for i in self.C}
        I[0] = np.round(max(self.status[0] + self.r - sum(self.q_[i, k] for k in self.K for i in self.C), 0), 4)
        U = {i: np.round(max(self.reveal[i] - (self.status[i] + sum(self.q_[i, k] for k in self.K)), 0), 4) for i in self.C}
        U[0] = 0

#        print(f"I = {I}")
#        print(f"U = {U}")

        self.fill_rate_hist[self.t + 1] = {'demand': sum(self.reveal[i] for i in self.C), 'lost_sales': sum(U[i] for i in self.C), 'fill_rate': 1-sum(U[i] for i in self.C)/sum(self.reveal[i] for i in self.C)}
        self.inventory_hist[self.t + 1] = I
        self.lost_sales_hist[self.t + 1] = U
        self.routes_hist[self.t + 1] = self.x_
        self.quantity_hist[self.t + 1] = self.q_
        self.status = status

#        print(f"status final= {status}")

        pass

    # siguiente periodo
    def next_time(self):
        self.t += 1
        pass

    # calcular costos totales
    def get_total_costs(self):
        self.inventory_costs = sum(self.h[i]*self.inventory_hist[t][i] for i in self.N for t in self.inventory_hist if t >= 6)
        self.lost_sales_costs = sum(self.a[i]*self.lost_sales_hist[t][i] for i in self.C for t in self.lost_sales_hist if t >= 6)
        self.routing_costs = 0

        for t in self.routes_hist:
            if t >= 6:
                for tup in self.routes_hist[t]:
                    self.routing_costs += self.routes_hist[t][tup]*self.c[tup[0], tup[1]]
                else:
                    pass

        self.total_costs = self.inventory_costs + self.lost_sales_costs + self.routing_costs
        pass


    def get_client_stats(self):
        def days(q):
            if q >= 0.01:
                return 1
            else:
                return 0

        self.inventory_node = {i: sum(self.h[i]*self.inventory_hist[t][i] for t in self.inventory_hist if t >= 6) for i in self.N}
        self.lost_sales_node = {i: sum(self.a[i]*self.lost_sales_hist[t][i] for t in self.inventory_hist if t >= 6) for i in self.C}
        self.quantity_node = {i: sum(self.quantity_hist[t][i, k] for k in self.K for t in self.quantity_hist if t >= 6) for i in self.C}
        self.visits_node = {i: sum(days(self.quantity_hist[t][i, k]) for k in self.K for t in self.quantity_hist if t >= 6) for i in self.C}
        self.lost_sales_days_node = {i: sum(days(self.lost_sales_hist[t][i]) for t in self.inventory_hist if t >= 6) for i in self.C}
        pass


#######################################################################################################
########################### PASO 6: DEFINIR DEFINIR OBJETO DE SIMULACION P2 ###########################
#######################################################################################################

class Simulation:
    def __init__(self):
        self.t = 0
        self.inventory = {}
        self.lost_sales = {}
        self.quantity = {}
        self.routing = {}


    def run(self, IRP):
        for t in range(3):
            IRP.action()
            IRP.reveal_demands()
            IRP.update_status()
            IRP.next_time()
        pass


##########################################################################################################
########################### PASO 7: DEFINIR FUNCIONES PARA LA TOMA DE DECISION ###########################
##########################################################################################################

# Modificacion 2

####ESTE ES EL IRP QUE HAY QUE RESOVER CON HEURISTICA!!!
def Solve_IRP_SS(sets, params, probs, status, t0 = 0, horizon = 3):

#    print("\nIRP EXACTO")
#    print(f"status = {status}")

    global T, C, N, A, K, x, y

    def IRPcut(model, where):

        def argmax(y, k, t, s):
            argmax = None
            maximum = -np.infty
            for i in s:
                if y[i, k, t] >= maximum:
                    argmax = i
                    maximum = y[i, k, t]
            return argmax

        def IRPsubtour(G):
            tour = nx.node_connected_component(G, 0)
            all_nodes = set(G.nodes())
            diff = all_nodes - tour
            return list(diff)

        if where == GRB.Callback.MIPSOL:

            x_ = model.cbGetSolution(x)
            y_ = model.cbGetSolution(y)

            for t in T:
                for k in K:
                    route_graph = nx.Graph()
                    route_graph.add_nodes_from([0] + [i for i in C if y_[i, k, t] >= 0.001])

                    if y_[0, k, t] >= 0.5:
                        for (i, j) in A:
                            if x_[i, j, k, t] >= 0.001:
                                route_graph.add_edge(i, j, capacity = x_[i, j, k, t])

                        subtour = IRPsubtour(route_graph)

                        if subtour:
                            for t1 in T:
                                for k1 in K:
                                    model.cbLazy(quicksum(x[i, j, k1, t1] for i in subtour for j in subtour if i != j) <= quicksum(y[i, k1, t1] for i in subtour) - y[argmax(y_, k, t, subtour), k1, t1])
                                    #model.cbLazy(quicksum(x[i, j, t] for i in subtour for j in subtour if i != j) <= len(subtour) - 1)

    T, C, N, A, K = [t0 + i + 1 for i in range(horizon)], sets['C'], sets['N'], sets['A'], sets['K']
    Q, cap, c, h, a, r = params['Q'], params['cap'], params['c'], params['h'], params['a'], params['r']
    safety_stock = params['ss']
    mu = probs['mu']
    t0 = t0
    T0 = [t0] + T

    m = Model()

    I = m.addVars(N, T0, name = 'I')
    q = m.addVars(C, K, T, name = 'q')
    x = m.addVars(A, K, T, vtype = GRB.BINARY, name = 'x')
    y = m.addVars(N, K, T, vtype = GRB.BINARY, name = 'y')
    u = m.addVars(C, T, name = 'u')

    obj = m.setObjective(quicksum(c[i, j]*x[i, j, k, t] for (i, j) in A for k in K for t in T) + quicksum(h[i]*I[i, t] for i in N for t in T) + quicksum(a[i]*u[i, t] for i in C for t in T))

#    print("restricciones de status")
#    print("r0  = m.addConstrs(I[i, t0] == status[i] for i in N)")
#    for i in N:
#        print(f"status[i] = {status[i]}")

    r0  = m.addConstrs(I[i, t0] == status[i] for i in N)
    r1  = m.addConstrs(I[i, t] == I[i, t - 1] + quicksum(q[i, k, t] for k in K) + u[i, t] - mu[i] for i in C for t in T)
    r2  = m.addConstrs(I[0, t] == I[0, t - 1] + r - quicksum(q[i, k, t] for i in C for k in K) for t in T)
    r3  = m.addConstrs(I[i, t - 1] + quicksum(q[i, k, t] for k in K) <= cap[i] for i in C for t in T)
    r4  = m.addConstrs(q[i, k, t] <= min(Q, cap[i])*y[i, k, t] for i in C for k in K for t in T)
    r5  = m.addConstrs(quicksum(q[i, k, t] for i in C) <= Q*y[0, k, t] for k in K for t in T)
    r6  = m.addConstrs(quicksum(x[i, j, k, t] for j in N if i != j) == y[i, k, t] for i in N for k in K for t in T)
    r7  = m.addConstrs(quicksum(x[j, i, k, t] for j in N if i != j) == y[i, k, t] for i in N for k in K for t in T)
    r8  = m.addConstrs(quicksum(y[i, k, t] for k in K) <= 1 for i in C for t in T)
    #r9 = subtour
    r10 = m.addConstrs(y[0, k, t] >= y[i, k, t] for i in C for k in K for t in T)
    r11 = m.addConstrs(x[i, j, k, t] <= y[i, k, t] for i in N for j in N for k in K for t in T if i != j)
    r12 = m.addConstrs(x[j, i, k, t] <= y[i, k, t] for i in N for j in N for k in K for t in T if i != j)
    r13 = m.addConstrs(y[0, k, t] <= y[0, k - 1, t] for k in K for t in T if k != 1)
    r14 = m.addConstrs(y[i, k, t] <= quicksum(y[j, k - 1, t] for j in C if j  < i) for i in C for k in K for t in T if k >= 2)
    r15 = m.addConstrs(I[i, t] >= min(cap[i] - mu[i], safety_stock[i, t-t0]) for i in C for t in T)

    m.params.OutputFlag = 0
    m.params.TimeLimit = 600
    m.params.LazyConstraints = 1
    m.update()

    m.optimize(IRPcut)

    x_ = {(i, j, k): round(x[i, j, k, t0 + 1].x) for (i, j) in A for k in K}
    q_ = {(i, k): q[(i, k, t0 + 1)].x for i in C for k in K}
    y_ = {(i, k): round(y[i, k, t0 + 1].x) for i in N for k in K}

    if horizon > 1:
        y_pred = {(i, k): round(y[i, k, t0 + 2].x, 2) for i in C for k in K}
    else:
        y_pred = {}

#    print(f"\nSolucion IRP Exacto :")
#    print(f"x_ : {x_}")
#    print(f"q_ : {q_}")
#    print(f"y_ : {y_}")
#    print(f"y_pred : {y_pred}")
    
    return x_, q_, y_, m.MIPGap, y_pred

############################################################## INICIO HEUR ##############################################################


def calculate_route_distance(route, cost_matrix):
    """
    Calculate the total cost of a given route.

    :param route: A list representing a route (sequence of nodes).
    :param cost_matrix: A dictionary representing the cost matrix with keys as node tuples.
    :return: Total cost of the route.
    """
    total_cost = 0
    for i in range(len(route) - 1):
        total_cost += cost_matrix.get((route[i], route[i + 1]), float('inf'))
    return total_cost



def report_instance(instance):

    print(f"\n### Instance Report ###")
    print(f"\nInstance Name: {instance['nombre']}")

    # Reporting horizon
    print(f"\nHorizon: {instance['horizon']}")

    # Reporting parameters
    print("\nParameters:")
    for param, value in instance['params'].items():
        print(f"  {param}: {value}")

    # Reporting sets
    print("\nSets:")
    for set_name, set_value in instance['sets'].items():
        print(f"  {set_name}: {set_value}")

    # Reporting demand
    print("\nDemand:")
    for customer, demand in instance['demanda'].items():
        print(f"  Customer {customer}: {demand}")

    # Reporting initial status
    print("\nInitial Status:")
    for node, status in instance['status'].items():
        node_label = "Depot" if node == 0 else f"Customer {node}"
        print(f"  {node_label}: {status}")

    # Optionally, report other elements as needed

    print(f"\n### Instance Report End ###")


def report_solution(solution):

    print("\n### Solution Report ###")
    print(f"\nName: {solution['nombre']}")

    print("\nRoutes:")
    for day, vehicles in solution['rutas'].items():
        print(f"  Day {day}:")
        for vehicle, route in vehicles.items():
            print(f"    Vehicle {vehicle}: {route}")

    print("\nDeliveries:")
    for customer, days in solution['deliveries'].items():
        print(f"  Customer {customer}: {days}")

    print("\nTotal Deliveries per Day:")
    for day, total in solution['deliveries_day'].items():
        print(f"  Day {day}: {total}")

    print("\nInventory Levels:")
    for node, inventory_levels in solution['inventarios'].items():
        node_label = "Depot" if node == 0 else f"Customer {node}"
        print(f"  {node_label}: {inventory_levels}")

    print("\nSafety Stock Deficit:")
    for customer, deficit in solution['deficit safety stock'].items():
        print(f"  Customer {customer}: {deficit}")

    print("\nRoute Lengths:")
    for day, lengths in solution['largos rutas'].items():
        print(f"  Day {day}:")
        for vehicle, length in lengths.items():
            print(f"    Vehicle {vehicle}: {length}")

    print(f"\nTotal Inventory Cost: {solution['costo inventario']}")
    print(f"Total Shortage Cost: {solution['costo deficit inventario']}")
    print(f"Total Vehicle Capacity Deficit Cost: {solution['costo deficit cap veh']}")
    print(f"Total Route Length: {solution['largo total']}")
    print(f"Total Cost: {solution['costo total']}")

def summary_solution(solution):
    print("\n### Solution Summary ###")
    print(f"\nName: {solution['nombre']}")

    print("\nRoutes:")
    for day, vehicles in solution['rutas'].items():
        print(f"  Day {day}:")
        for vehicle, route in vehicles.items():
            print(f"    Vehicle {vehicle}: {route}")

    print(f"\nTotal Inventory Cost: {solution['costo inventario']}")
    print(f"Total Shortage Cost: {solution['costo deficit inventario']}")
    print(f"Total Vehicle Capacity Deficit Cost: {solution['costo deficit cap veh']}")
    print(f"Total Route Length: {solution['largo total']}")
    print(f"Total Cost: {solution['costo total']}")
    print("\n### Solution Summary End ###")

def construye_solucion(instancia, routes):
    
    #Construye un diccionario que contiene las rutas, inventarios, y costos de una solución

    solucion = {}
    solucion['nombre'] = "NN"
    solucion['rutas'] = routes

    update_solution_from_routes(solucion, instancia)

    return solucion

def update_solution_from_routes(solucion, instancia):

    horizon = instancia['horizon']
    T, C, N, A, K = [i + 1 for i in range(horizon)], instancia['sets']['C'], instancia['sets']['N'], instancia['sets']['A'], instancia['sets']['K']
    Q, cap, c, h, a, r, ss = instancia['params']['Q'], instancia['params']['cap'], instancia['params']['c'], instancia['params']['h'], instancia['params']['a'], instancia['params']['r'], instancia['params']['ss']
    mean_demand = instancia['demanda']
    status = instancia['status']
    routes = solucion['rutas']

#    print(f"routes para update : {routes}")

    #Calcula deliveries, inventarios y déficits

    delivery = {cliente: {day: 0 for day in [1, 2 ,3]} for cliente in C}  ### HARD CODED, OJO!
    inventory = {cliente: {day: 0 for day in [0, 1, 2 ,3]} for cliente in N}
    deficit = {cliente: {day: 0 for day in [1, 2 ,3]} for cliente in C}
    delivery_day = {day: 0 for day in [1, 2 ,3]}

    for day in T:
        for veh in K:
            for cliente, cantidad in routes[day][veh]:
                delivery[cliente][day] = delivery[cliente][day] + cantidad
                delivery_day[day] = delivery_day[day] + cantidad
    
    solucion['deliveries'] = delivery
    solucion['deliveries_day'] = delivery_day

    for cliente in N:

        inventory[cliente][0] = status[cliente]

        for day in T:
            if cliente != 0:
                inventory[cliente][day] = max(inventory[cliente][day-1] - mean_demand[cliente] + delivery[cliente][day], 0)
                if inventory[cliente][day] < ss[(cliente, day)]:
                    deficit[cliente][day] = mean_demand[cliente] - max(inventory[cliente][day] - ss[(cliente, day)],0)
            else:
                inventory[cliente][day] = max(inventory[cliente][day-1] - delivery_day[day] + r, 0)

    solucion['inventarios'] = inventory
    solucion['deficit safety stock'] = deficit

    for day in T:
        inventory[cliente][day] = max(inventory[cliente][day-1] - mean_demand[cliente] + delivery[cliente][day], 0)
        if inventory[cliente][day] < ss[(cliente, day)]:
            deficit[cliente][day] = mean_demand[cliente] - max(inventory[cliente][day] - ss[(cliente, day)],0)

    #Costos de rutas

    largos_rutas = {day: {veh: 0 for veh in K} for day in T}

    for day in T:
        for veh in K:
            ruta_diaria = routes[day][veh]
            dist_ruta = 0

            #Si la ruta tiene 1 cliente o más, considerar dists desde y hacia el depot
            if len(ruta_diaria) >= 1:
                dist_ruta = dist_ruta + c[(0, ruta_diaria[0][0])] + c[(ruta_diaria[-1][0], 0)]

            #Si la ruta tiene 2 clientes o más, considerar además distancias entre clientes
            if len(ruta_diaria) >= 2:
                for pos in range(1, len(ruta_diaria)):
                    dist_ruta = dist_ruta + c[(ruta_diaria[pos - 1][0], ruta_diaria[pos][0])]

            largos_rutas[day][veh] = dist_ruta

    solucion['largos rutas'] = largos_rutas

    #Costos de inventario

    inv_cost = 0
    def_cost = 0

    for day in T:
        for node in N:
            inv_cost = inv_cost + h[node] * inventory[node][day]
        for cliente in C:
            def_cost = def_cost + 10 * deficit[cliente][day] #### OJO QUE AQUI ESTA LA PENALIDAD POR DEFICIT
    
    solucion['costo inventario'] = inv_cost
    solucion['costo deficit inventario'] = def_cost

    #Penalidad por deficit de capacidad vehículos
            
    total_deliveries = {}
    vehicle_cap_deficit = {}

    for day, vehicles in routes.items():
        total_deliveries[day] = {}
        vehicle_cap_deficit[day] = {}
        for vehicle, deliveries in vehicles.items():
            total_quantity = sum(quantity for _, quantity in deliveries)
            total_deliveries[day][vehicle] = total_quantity
            vehicle_cap_deficit[day][vehicle] = max(total_deliveries[day][vehicle] - Q, 0)

    solucion['deficit cap veh'] = vehicle_cap_deficit

    veh_deficit_cost = 0
    largo_total = 0

    for day in T:
        for vehicle in K:
            veh_deficit_cost = veh_deficit_cost + 1000 * vehicle_cap_deficit[day][vehicle] #### OJO QUE AQUI ESTA LA PENALIDAD POR DEFICIT
            largo_total = largo_total + largos_rutas[day][vehicle]

    solucion['costo deficit cap veh'] = veh_deficit_cost
    solucion['largo total'] = largo_total
    solucion['costo total'] = veh_deficit_cost + def_cost + inv_cost + largo_total

    # Variables para el período 1

    # Initialize entries
    solucion['x_'] = {(i, j, k): 0 for (i, j) in A for k in K}
    solucion['q_'] = {(i, k): 0 for i in C for k in K}
    solucion['y_'] = {(i, k): 0 for i in N for k in K}
    solucion['y_pred'] = {(i, k): 0 for i in C for k in K}

    # Update based on routes for day 1
    day = 1

    for vehicle, route in solucion['rutas'][day].items():
        prev_node = 0  # Starting from the depot
        for customer, quantity in route:
            # Update x_
            solucion['x_'][(prev_node, customer, vehicle)] = 1

            # Update q_ and y_
            if customer != 0:  # If not the depot
                solucion['q_'][(customer, vehicle)] = quantity
                solucion['y_'][(customer, vehicle)] = 1

            prev_node = customer

        # Include the return arc to the depot
        if route:
            solucion['x_'][(route[-1][0], 0, vehicle)] = 1  # Last customer to depot

    for vehicle in K:
        if (solucion['rutas'][1][vehicle]):
            solucion['y_'][(0, vehicle)] = 1

    # Update based on routes for day 2

    day = 2

    for vehicle, route in solucion['rutas'][day].items():
        prev_node = 0  # Starting from the depot
        for customer, quantity in route:
            if customer != 0:  # If not the depot
                solucion['y_pred'][(customer, vehicle)] = 1
            prev_node = customer

    return solucion

def plan_visit_dates(instancia):

    horizon = instancia['horizon']
    T, C, N, A, K = [i + 1 for i in range(horizon)], instancia['sets']['C'], instancia['sets']['N'], instancia['sets']['A'], instancia['sets']['K']
    Q, cap, c, h, a, r, ss = instancia['params']['Q'], instancia['params']['cap'], instancia['params']['c'], instancia['params']['h'], instancia['params']['a'], instancia['params']['r'], instancia['params']['ss']
    mean_demand = instancia['demanda']
    initial_inventory = instancia['status']

    visit_plan = {customer: [] for customer in C}

    for customer in C:

        inv = initial_inventory[customer]

        for day in range(1, horizon + 1):
            if (inv - mean_demand[customer]) <= ss[(customer, day)]:
                inv = cap[customer] - mean_demand[customer]
                visit_plan[customer].append(day)
            else:
                inv = inv - mean_demand[customer]

    return visit_plan

def compute_required_deliveries(instancia, visit_plan):

    horizon = instancia['horizon']
    T, C, N, A, K = [i + 1 for i in range(horizon)], instancia['sets']['C'], instancia['sets']['N'], instancia['sets']['A'], instancia['sets']['K']
    Q, cap, c, h, a, r, ss = instancia['params']['Q'], instancia['params']['cap'], instancia['params']['c'], instancia['params']['h'], instancia['params']['a'], instancia['params']['r'], instancia['params']['ss']
    mean_demand = instancia['demanda']
    initial_inventory = instancia['status']

    # Initialize inventory levels
    daily_inventories = {customer: [initial_inventory[customer]] for customer in C}
    daily_deliveries = {customer: [0] for customer in C}

    # Simulate each day
    for day in range(1, horizon + 1):  # Day 0 is the initial inventory

        for customer in C:

            # Check if there is a visit to the customer on current day:
            if day in visit_plan[customer]:
                daily_inventories[customer].append(cap[customer] - mean_demand[customer])
                daily_deliveries[customer].append(cap[customer] - daily_inventories[customer][day-1])
            else:
                daily_inventories[customer].append(daily_inventories[customer][-1] - mean_demand[customer])
                daily_deliveries[customer].append(0)

    return daily_deliveries

def calculate_total_daily_deliveries(daily_deliveries_amount):

    total_daily_deliveries = {}

    for customer, deliveries in daily_deliveries_amount.items():
        for day, amount in enumerate(deliveries):
            if day in total_daily_deliveries:
                total_daily_deliveries[day] += amount
            else:
                total_daily_deliveries[day] = amount

    return total_daily_deliveries

def build_daily_vrp_solutions(instancia, daily_deliveries_amount):

    horizon = instancia['horizon']
    T, C, N, A, K = [i + 1 for i in range(horizon)], instancia['sets']['C'], instancia['sets']['N'], instancia['sets']['A'], instancia['sets']['K']
    Q, cap, c, h, a, r, ss = instancia['params']['Q'], instancia['params']['cap'], instancia['params']['c'], instancia['params']['h'], instancia['params']['a'], instancia['params']['r'], instancia['params']['ss']
    mean_demand = instancia['demanda']
    initial_inventory = instancia['status']

    vrp_horizon = range(1, 4)

    vehicle_capacity = Q
    vehicles = K
    daily_routes = {day: {vehicle: [] for vehicle in vehicles} for day in vrp_horizon}

    for day in vrp_horizon:
        customers_to_visit = [customer for customer, deliveries in daily_deliveries_amount.items() if deliveries[day] > 0]
        demands = {customer: daily_deliveries_amount[customer][day] for customer in customers_to_visit}
        remaining_capacities = {vehicle: vehicle_capacity for vehicle in vehicles}

        for vehicle in vehicles:
            for customer in list(customers_to_visit):
                if demands[customer] <= remaining_capacities[vehicle]:
                    # Append a tuple of (customer, quantity)
                    daily_routes[day][vehicle].append((customer, demands[customer]))
                    remaining_capacities[vehicle] -= demands[customer]
                    customers_to_visit.remove(customer)

        # Assign remaining customers to the vehicle with the most available capacity
        for customer in customers_to_visit:
            vehicle_with_most_space = max(remaining_capacities, key=remaining_capacities.get)
            if demands[customer] <= remaining_capacities[vehicle_with_most_space]:
                daily_routes[day][vehicle_with_most_space].append((customer, demands[customer]))
                remaining_capacities[vehicle_with_most_space] -= demands[customer]

    return daily_routes



####################################################################################################
########################### PASO 8: CORRER LA SIMULACION #################3#########################
####################################################################################################

# Modificacion 3

#ALGORITMOS QUE SE VAN A CORRER
option_names = {1: 'IRP(1)', 2: 'IRP(3)', 3: 'IRP(6)',
                4: 'IRP(1)-SS', 5: 'IRP(3)-SS', 6: 'IRP(6)-SS',
                7: 'IRP(1)-RP(1)', 8: 'IRP(1)-RP(2)', 9: 'IRP(1)-RP(3)', 10: 'IRP(1)-RP(4)',
                11: 'IRP(2)-RP(1)', 12: 'IRP(2)-RP(2)', 13: 'IRP(2)-RP(3)', 14: 'IRP(2)-RP(4)',
                15: 'IRP(3)-RP(1)', 16: 'IRP(3)-RP(2)', 17: 'IRP(3)-RP(3)', 18: 'IRP(3)-RP(4)',
                19: 'IRP(1)-RP(1)-SS', 20: 'IRP(1)-RP(2)-SS', 21: 'IRP(1)-RP(3)-SS', 22: 'IRP(1)-RP(4)-SS',
                23: 'IRP(2)-RP(1)-SS', 24: 'IRP(2)-RP(2)-SS', 25: 'IRP(2)-RP(3)-SS', 26: 'IRP(2)-RP(4)-SS',
                27: 'IRP(3)-RP(1)-SS', 28: 'IRP(3)-RP(2)-SS', 29: 'IRP(3)-RP(3)-SS', 30: 'IRP(3)-RP(4)-SS', 31: 'HEUR'}

random.seed('EXP1')

cols = ['ID', 'CID', '|C|', '|K|', 'OPTION', 'INVENTORY', 'LOST SALES', 'ROUTING', 'TOTAL COST', 'TIME']
results = {col: [] for col in cols}

dis_cols = ['ID', 'CID', '|C|', '|K|', 'OPTION', 'INVENTORY', 'LOST SALES', 'ROUTING', 'TOTAL COST', 'QUANTITY', 'VISITS', 'PREDICTIONS','FULL DEMAND', 'LOST DEMAND', 'FILL RATE', 'MIP GAP','PERIOD']
dis_results = {col: [] for col in dis_cols}
def simular_ejecucion_P_deterministico(grafo_inicial, dem_historico, capacidad, tipo_demanda = dem, T = 365, d=30, F=7)
    
    G0 = grafo_inicial.copy()
    matriz_dst = calcular_matriz_dist(G0)
    ubicaciones = list(G0.nodes()) # Lista de ubicaciones
    inventarios = [[G0.nodes(data=True)[i]['Inv'] for i in ubicaciones]] # Lista de inventarios
    h = [G0.nodes(data=True)[i]['h'] for i in ubicaciones] # Lista de costos de inventario
    rutas = {t : [] for t in range(T)} # Lista de rutas

<<<<<<< HEAD
    
=======
    inventario_total = []
    perdidas = []
    c_rutas =[]
    demandas_efectivas = []
    costo_rutas = 0
    costo_SO = 0

    #self.T = unpickled['Sets']['T'] #periodos
    #self.C = unpickled['Sets']['C'] #clientes
    #self.N = unpickled['Sets']['N'] #nodos
    #self.A = unpickled['Sets']['A'] #arcos
    #self.K = unpickled['Sets']['K'] #vehiculos

    #Parameters
    #h = unpickled['Param']['h'] #costo de inventario
    #self.Q = unpickled['Param']['Q'] #capacidad de vehiculo
    #self.cap = unpickled['Param']['cap'] #capacidad de bodega
    #self.pos = unpickled['Pos']['pos'] #matriz de posiciones en el plano cartesiano
    #self.c = unpickled['Pos']['cost'] #matriz de costos de transporte
    #self.h = {i: h[i] for i in self.N for t in self.T} #costo de inventario
    #self.a = unpickled['Param']['a'] #costo de demanda insatisfecha
    #self.r = unpickled['Param']['r'] #cantidad de unidades que llegan al depot en cada periodo
    for t in range(T):
        print('\n')
        mu_demanda = [np.mean(dem_historico[nodo]) for nodo in dem_historico.keys()]    
        sd_demanda = [np.std(dem_historico[nodo]) for nodo in dem_historico.keys()]
        pronostico = {int(nodo[2:]): pronostico_SEDA(
                                    dem_historico[nodo], T = F, pron = True, alpha=0.2, beta=0.1, theta=0.5)[0]
                                    for nodo in dem_historico.keys()}
        # print(pronostico)
        pronostico = adaptar_pron(pronostico, F)

        # ahora incluimos el modelo de optimización
        params = {}
        params['Q'] = capacidad
        params['cap'] = {i: G.nodes(data=True)[i]['Up'] for i in ubicaciones}
        params['c'] = matriz_dst
        params['h'] = h
        params['a'] = 10 # REVISAR ESTE COSTO
        params['r'] = 30 # revisar esto también, HAY QUE CAMBIARLO
        params['pos'] = {i: G.nodes(data=True)[i]['pos'] for i in ubicaciones}

        sets = {}
        sets['T'] = [t + t_i for t_i in range(F)]
        sets['C'] = [i for i in ubicaciones if i != 'N_0']
        sets['N'] = ubicaciones
        sets['A'] = [e for e in G0.edges()]
        sets['K'] = 1

        status = {i: G0.nodes(data=True)[i]['Inv'] for i in ubicaciones} # MODIFICAR!!!

        PROBS = {}

        solucion = Solve_IRP_SS(sets, params, probs, status, t0 = t, horizon = F)

        # print(f"solucion = {solucion}")

        cr = calcular_largo_ruta(ruta_P, distancias)
        c_rutas.append(cr)
        costo_rutas += cr

        if ruta_P != [] and ruta_P != None and ruta_P != ['N_0']:
             ruta_P += ['N_0']
             G0, stock = ejecutar_ruta(G0, ruta_P, distancias)
        
        elif ruta_P == ['N_0'] or ruta_P == ['N_0','N_0']:
            ruta_P = []

        rutas[t] = ruta_P
        # print(f"Ruta {t}: ", ruta_P)
        # visitas_proactiva = proactiva_inventario(G0, tolerancia = 0.2, dist = 'n', mu = 0, sigma = 0.1, M = 1000)

        G0, demanda, insatisfecho = realizacion_demanda_modificada(G0, dist = tipo_demanda, T=T, demandas_in=dem_historico, d=d, t=t)
        demandas_efectivas.append(demanda)
        costo_SO += insatisfecho*1
        d_total += sum(demanda.values())
        inventarios = [G0.nodes(data=True)[i]['Inv'] for i in ubicaciones if i != 'N_0']
        inventario_total.append(sum(inventarios))
        perdidas.append(insatisfecho)

        print(f'Tiempo: {t} | Ruta: {ruta_P} | costo_SO: {insatisfecho*1} | costo_r: {cr}')
        #Actualizo demandas
        for nodo in ubicaciones:
            if nodo != 'N_0':
                dem_historico[nodo].append(demanda[nodo]) 
                # HAY QUE VER CÓMO SE COMPORTA ESTO CON LA DEMANDA MODIFICADA

    # print('\n')
    # print("Inventario final: ")
    # for nodo in G0.nodes(data=True):
    #     print(nodo[0],nodo[1]['Inv'])
    print(f'F = {F}, Demanda perdida total: {sum(perdidas)} | Demanda perdida promedio: {sum(perdidas)/T}')        
    print(f'Costo ruta: {costo_rutas} | Costo SO: {costo_SO}')
    # graficar_rutas(rutas, G0)
    costos = (perdidas, c_rutas)
    return rutas, perdidas, inventario_total, costos, demandas_efectivas
>>>>>>> 7958b871f0e863ac9f6a46fd3d57cc9f2c2e5468
    
    






<<<<<<< HEAD

=======
"""
>>>>>>> 7958b871f0e863ac9f6a46fd3d57cc9f2c2e5468

    for n_clients in [10]:
        for n_vehicles in [2]:
            op = 5

            for vid in [0]:

                dir = 'Archivos pickle/'
                name = 'I{}-K{}-V{}'.format(n_clients, n_vehicles, vid + 1)
                rd = pd.read_pickle(dir + 'SDEMAND-{}.pickle'.format(name))
                cost_params = pd.read_pickle(dir + 'cost_params_ocf-I{}-K{}-V{}.pickle'.format(n_clients, n_vehicles, vid + 1))

                #Simulacion
                problem = IRP(dir + name + '.pickle', rd, op, cost_params)
                simulation = Simulation()
                toc = time.time()
                simulation.run(problem)
                tic = time.time()

                problem.get_total_costs()

                #Exportar resultados
                results['ID'].append(vid + 1)
                results['CID'].append(cid + 1)
                results['|C|'].append(len(problem.C))
                results['|K|'].append(len(problem.K))
                results['OPTION'].append(option_names[op])
                results['INVENTORY'].append(problem.inventory_costs)
                results['LOST SALES'].append(problem.lost_sales_costs)
                results['ROUTING'].append(problem.routing_costs)
                results['TOTAL COST'].append(problem.total_costs)
                results['TIME'].append(tic-toc)
                df = pd.DataFrame(data=results)
                df.to_csv('Results/d_reults_I{}.csv'.format(n_clients), decimal = ',', sep = '\t')



                df = pd.DataFrame(data=dis_results)
                df.to_csv('Results/d_dis_reults_I{}.csv'.format(n_clients), decimal = ',', sep = '\t')
<<<<<<< HEAD

=======
"""
>>>>>>> 7958b871f0e863ac9f6a46fd3d57cc9f2c2e5468
# %%
