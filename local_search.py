'''
FUNCIONES PARA LA IMPLEMENTACIÓN DE LA POLÍTICA PROACTIVA CON LOCAL SEARCH
'''


import numpy as np
from numpy import random as rd
import time
from solvers_listos import *
from funciones import *



def eliminar_duplicados(ruta):
    unique_nodes = []
    [unique_nodes.append(node) for node in ruta if node not in unique_nodes]
    return unique_nodes

def random_reverse(rutas, nodos = [] ,iters = 1):
    for _ in range(iters):
        id = rd.choice(len(rutas))
        ruta = rutas[id]
        if len(ruta) > 2:
            i, j = rd.choice(range(1, len(ruta)), 2, replace=False)
            ruta[i], ruta[j] = ruta[j], ruta[i]
            # print(f'{id}: {ruta}')
        rutas[id] = ruta
    return rutas

def random_swap(rutas, nodos = [], iters=1):
    for _ in range(iters):
        id1,id2 = rd.choice(len(rutas), 2, replace = False)
        ruta1, ruta2 = rutas[id1], rutas[id2]
        # print(f'{id1}: {ruta1} | {id2}: {ruta2}')
        if len(ruta1) > 3 and len(ruta2)>3:
            bloque1 = rd.choice(range(1, len(ruta1)-1))  
            par1 = (ruta1[bloque1], ruta1[bloque1 + 1])
            bloque2 = rd.choice(range(1, len(ruta2)-1))
            par2 = (ruta2[bloque2], ruta2[bloque2 + 1])
            # print(f'R{id1}-Bloque1: {par1} x R{id2}-Bloque2: {par2}')
            ruta1[bloque1], ruta1[bloque1+1] = par2[0], par2[1] 
            ruta2[bloque2], ruta2[bloque2+1] = par1[0], par1[1]
            ruta1, ruta2 = eliminar_duplicados(ruta1), eliminar_duplicados(ruta2)
        rutas[id1] = ruta1
        rutas[id2] = ruta2
    return rutas

def random_move(rutas, nodos = [], iters = 1):
    for _ in range(iters):
        id = rd.choice(len(rutas))
        ruta = rutas[id]
        if len(ruta) > 3:
            i, j = rd.choice(range(1, len(ruta)), 2, replace=False)
            node = ruta.pop(i)
            # print(f'Muevo el nodo {node} de la ruta {id} de la posición {i} a la posición {j}')
            ruta.insert(j, node)
        rutas[id] = ruta
    return rutas

def random_insert(rutas, nodos, iters = 1):
    for _ in range(iters):
        id = rd.choice(len(rutas))
        ruta = rutas[id]
        nodo = rd.choice(nodos)
        if ruta != [] and len(ruta)>2:
            pos = rd.choice(range(1, len(ruta)))
            # print(f'Inserto el nodo {nodo} en la ruta {id} en la posición {pos}')
            ruta.insert(pos, nodo)
            ruta = eliminar_duplicados(ruta)
            rutas[id] = ruta
    return rutas

def random_remove(rutas, nodos = [], iters = 1):
    for _ in range(iters):
        id = rd.choice(len(rutas))
        ruta = rutas[id]
        if ruta != [] and len(ruta) > 1:
            pos = rd.choice(range(1, len(ruta)))
            # print(f'Elimino el nodo {ruta[pos]} de la ruta {id} en la posición {pos}')
            ruta.pop(pos)
            rutas[id] = ruta
    return rutas

def operacion_random(rutas, nodos):
    ops = [random_insert, random_remove, random_move, random_swap, random_reverse]
    op = rd.choice(ops)
    # print(f'Operación: {op.__name__}')
    return op(rutas,nodos)    

def costo_total(rutas, distancias):
    costo = 0
    for ruta in rutas.values():
        if ruta != []:
            for i in range(len(ruta) - 1):
                costo += distancias[ruta[i]][ruta[i + 1]]
        else:
            costo += 0
    return costo

def Local_Search(G, ruta_0, demandas, distancias, cap, F, n_restarts = 10, n_iters = 5):
    
    # time_limit = 10
    # t0 = time.time()
    best = ruta_0.copy()
    # print(best)
    best_eval = costo_total(best, distancias)
    nodos = G.nodes()

    # while time.time() - t0 < time_limit:
    for t in range(n_restarts):
        new = best.copy()
        for k in range(n_iters):
            new = operacion_random(new, nodos)

        new_eval = costo_total(new, distancias)
        if new_eval < best_eval:
            best, best_eval = new, new_eval
            # print(f'Restart {t}, best: {best_eval}')
    # print(f'Best: {best}')
    return best, best_eval

def realizacion_demanda_LS(G, demandas, ruido = 0.05):
    """
    Función que simula la demanda de los locales para un determinado periodo.
    """
    grafo = G.copy()
    insatisfecho = 0
    for nodo in grafo.nodes(data=True):
        if nodo[0] != 'N_0':
            dem = demandas[int(nodo[0][2:])]
            if dem <= grafo.nodes[nodo[0]]['Inv']:
                grafo.nodes[nodo[0]]['Inv'] -= dem

            else:
                grafo.nodes[nodo[0]]['Inv'] = 0
                insatisfecho += dem - grafo.nodes[nodo[0]]['Inv']
    # print(dems)
    # for nodo in grafo.nodes(data=True):
        # print(nodo[0],nodo[1]['Inv'])

    return grafo, insatisfecho

def adaptar_pron(prono, F):
    dict_pro = {}
    for t in range(F):
        dict_pro[t] = {nodo: prono[nodo]  for nodo in prono.keys()}
    return dict_pro

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

def ruteo_LS(H, distancias, pron_demandas, cap, F, mu, sd):
    grafo = H.copy()

    rutas = {}
    for t in range(F):
        '''
        Resolver el problema de ruteo para el periodo t
        Demanda = pronostico_demandas[t]
        '''
        demandas_t = pron_demandas[t]
        visitas_NN = reaccion_inventario(grafo, mu, sd)
        # print(visitas_NN)
        if sum(visitas_NN.values()) == 0:
            # print("No hay locales que visitar")
            rutas[t] = []

        else:    
            ruta_R = nearest_neighbor(grafo, distancias, disponibilidad=visitas_NN)[-1] #devuelve la ruta a realizar
            rutas[t] = ruta_R
            # print(f"Ruta {t}: ", ruta_R )
            grafo, stock = ejecutar_ruta(grafo, ruta_R, distancias)   
        grafo, insatisfecho = realizacion_demanda_LS(grafo, demandas_t)

    rutas, costo = Local_Search(grafo, rutas, demandas_t, distancias, cap, F)

    ruta_LS = eliminar_duplicados(rutas[0])

    return ruta_LS

def simular_ejecucion_P_LS(grafo_inicial, cap, dem_historico, T=1, F=1):
    # Inicializar variables     
    # ---------------------
    G0 = grafo_inicial.copy()
    distancias = calcular_matriz_dist(G0)
    ubicaciones = list(G0.nodes()) # Lista de ubicaciones
    inventarios = [G0.nodes(data=True)[i]['Inv'] for i in ubicaciones] # Lista de inventarios
    h = [G0.nodes(data=True)[i]['h'] for i in ubicaciones] # Lista de costos de inventario
    # dem_historico = simular_demanda_previa(G0, dist = 'n', T=1000) 
    d_total = 0
    rutas = {t : None for t in range(T)} # Lista de rutas
    inventario_total = []
    perdidas = []

    # print("Inventario inicial: ")
    # for nodo in G0.nodes(data=True):
    #         print(nodo[0],nodo[1]['Inv'])
    # print("\n")

    for t in range(T):
        # print('\n')
        mu_demanda = [np.mean(dem_historico[nodo]) for nodo in dem_historico.keys()]    
        sd_demanda = [np.std(dem_historico[nodo]) for nodo in dem_historico.keys()]
        pronostico = {int(nodo[2:]): pronostico_SEDA(
                                    dem_historico[nodo], T = F, pron = True, alpha=0.2, beta=0.1, theta=0.5)[0]
                                    for nodo in dem_historico.keys()}
        # print(pronostico)
        pronostico = adaptar_pron(pronostico, F)

        ruta_P = ruteo_LS(H = G0, distancias = distancias, pron_demandas = pronostico,
                           cap = cap, F = F, mu = mu_demanda, sd =sd_demanda)
        

        if ruta_P != [] and ruta_P != None and ruta_P != ['N_0']:
             ruta_P += ['N_0']
             G0, stock = ejecutar_ruta(G0, ruta_P, distancias)
        
        elif ruta_P == ['N_0']:
            ruta_P = []

        rutas[t] = ruta_P
        # print(f"Ruta {t}: ", ruta_P)
        # visitas_proactiva = proactiva_inventario(G0, tolerancia = 0.2, dist = 'n', mu = 0, sigma = 0.1, M = 1000)

        
        G0, demanda, insatisfecho = realizacion_demanda(G0)
        d_total += demanda
        inventarios = [G0.nodes(data=True)[i]['Inv'] for i in ubicaciones if i != 'N_0']
        inventario_total.append(sum(inventarios))
        perdidas.append(insatisfecho)

        #Actualizo demandas
        for nodo in ubicaciones:
            if nodo != 'N_0':
                dem_historico[nodo].append(demanda[nodo])


    # print('\n')
    # print("Inventario final: ")
    # for nodo in G0.nodes(data=True):
    #     print(nodo[0],nodo[1]['Inv'])
    print(f'F = {F}, Demanda perdida total: {sum(perdidas)} | Demanda perdida promedio: {sum(perdidas)/T}')        

    # graficar_rutas(rutas, G0)
    return rutas, perdidas, inventario_total

# rutas, perdidas, inventarios = simular_ejecucion_P_LS(grafo_inicial = G, T = 365, F = 10, cap = 871)