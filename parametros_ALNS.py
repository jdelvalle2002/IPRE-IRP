G, ubis, cap_tpte, info_locales = crear_grafo_inicial(archivo= 'IRP1.xlsx' ,plot=False)
matriz_dst = calcular_matriz_dist_alns(G)
demandas = simular_demanda_diaria(list(info_locales["r"]), dist="n")
print(demandas)

degree_of_destruction = 0.2
nodes_to_destroy = int(np.ceil(degree_of_destruction * len(G.nodes)))