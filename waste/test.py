# iterated local search of the ackley objective function
from numpy import asarray
from numpy import exp
from numpy import sqrt
from numpy import cos
from numpy import e
from numpy import pi
from numpy.random import randn
from numpy.random import rand
from numpy.random import seed
import random
import time

# objective function
def objective(v):
	x, y = v
	return -20.0 * exp(-0.2 * sqrt(0.5 * (x**2 + y**2))) - exp(0.5 * (cos(2 * pi * x) + cos(2 * pi * y))) + e + 20

# check if a point is within the bounds of the search
def in_bounds(point, bounds):
	# enumerate all dimensions of the point
	for d in range(len(bounds)):
		# check if out of bounds for this dimension
		if point[d] < bounds[d, 0] or point[d] > bounds[d, 1]:
			return False
	return True

# hill climbing local search algorithm
def hillclimbing(objective, bounds, n_iterations, step_size, start_pt):
	# store the initial point
	solution = start_pt
	# evaluate the initial point
	solution_eval = objective(solution)
	# run the hill climb
	for i in range(n_iterations):
		# take a step
		candidate = None
		while candidate is None or not in_bounds(candidate, bounds):
			candidate = solution + randn(len(bounds)) * step_size
		# evaluate candidate point
		candidte_eval = objective(candidate)
		# check if we should keep the new point
		if candidte_eval <= solution_eval:
			# store the new point
			solution, solution_eval = candidate, candidte_eval
	return [solution, solution_eval]

# iterated local search algorithm
def iterated_local_search(objective, bounds, n_iter, step_size, n_restarts, p_size):
	# define starting point
	best = None
	while best is None or not in_bounds(best, bounds):
		best = bounds[:, 0] + rand(len(bounds)) * (bounds[:, 1] - bounds[:, 0])
	# evaluate current best point
	best_eval = objective(best)
	# enumerate restarts
	for n in range(n_restarts):
		# generate an initial point as a perturbed version of the last best
		start_pt = None
		while start_pt is None or not in_bounds(start_pt, bounds):
			start_pt = best + randn(len(bounds)) * p_size
		# perform a stochastic hill climbing search
		solution, solution_eval = hillclimbing(objective, bounds, n_iter, step_size, start_pt)
		# check for new best
		if solution_eval < best_eval:
			best, best_eval = solution, solution_eval
			print('Restart %d, best: f(%s) = %.5f' % (n, best, best_eval))
	return [best, best_eval]

# seed the pseudorandom number generator
seed(1)
# define range for input
bounds = asarray([[-5.0, 5.0], [-5.0, 5.0]])
# define the total iterations
n_iter = 1000
# define the maximum step size
s_size = 0.05
# total number of random restarts
n_restarts = 30
# perturbation step size
p_size = 1.0
# perform the hill climbing search
best, score = iterated_local_search(objective, bounds, n_iter, s_size, n_restarts, p_size)
print('Done!')
print('f(%s) = %f' % (best, score))









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
    dist = distancias[pred][node] + distancias[node][succ] - distancias[pred][succ]
    return dist

def hillclimbing(objective, bounds, n_iterations, step_size, start_pt):
	# store the initial point
	solution = start_pt
	# evaluate the initial point
	solution_eval = objective(solution)
	# run the hill climb
	for i in range(n_iterations):
		# take a step
		candidate = None
		while candidate is None or not in_bounds(candidate, bounds):
			candidate = solution + randn(len(bounds)) * step_size
		# evaluate candidate point
		candidte_eval = objective(candidate)
		# check if we should keep the new point
		if candidte_eval <= solution_eval:
			# store the new point
			solution, solution_eval = candidate, candidte_eval
	return [solution, solution_eval]


def random_reverse(rutas):
    for _ in range(2):
        route = random.choice(rutas)
        if len(route) > 2:
            i, j = random_state.choice(range(1, len(route) - 1), 2, replace=False)
            route[i], route[j] = route[j], route[i]
    return swapped

rutas = [[1,5,3,4],[2,6,1,7,5], [2,6,1,7,5]]


def random_insertion(rutas):
    for _ in range(2):
        route = random_state.choice(rutas)
        if len(route) > 2:
            i, j = random_state.choice(range(1, len(route) - 1), 2, replace=False)
            node = route.pop(i)
            route.insert(j, node)
    return swapped



def Local_Search(G, ruta_0, demandas, distancias, cap, F, n_restarts = 10):
    best = ruta_0
    time_limit = 10
    t0 = time.time()
    best_eval = objective(best)
    while time.time() - t0 < time_limit:
        
        # generate an initial point as a perturbed version of the last best
        start_pt = None
        while start_pt is None or not in_bounds(start_pt, bounds):
            start_pt = best + randn(len(bounds)) * p_size
        # perform a stochastic hill climbing search
        solution, solution_eval = hillclimbing(objective, bounds, n_iter, step_size, start_pt)
        # check for new best
        if solution_eval < best_eval:
            best, best_eval = solution, solution_eval
            print('Restart %d, best: f(%s) = %.5f' % (n, best, best_eval))

    return [best, best_eval]