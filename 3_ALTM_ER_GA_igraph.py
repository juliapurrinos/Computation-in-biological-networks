import numpy as np
from numpy.random import choice
import networkx as nx
import random
import igraph as ig
import os
import sys, argparse


directory_path = "C:/Users/ba03/OneDrive/Desktop/Programming/FYP/output/z/2/"


'''
activation iteration model: 
iterates through each possible network and for each possible network it iterates through 
each possible input pattern then saves the data in a numpy array.
Each for loop of network generation = 1 trial
'''

# def Correct_Usage():
#         print("Correct Way to run program is:")
#         print("progname -n number of nodes -t number of trials -z zvalue -a antagonism -o outputname")
#         sys.exit(1)


# # ==== Main program starts here

# parser = argparse.ArgumentParser()
# parser.add_argument("-n", "--nodes"    , dest = "nodes",   default = "none", help="Nodes")
# parser.add_argument("-t", "--trials"  , dest = "trials", default = "none", help="Trials")
# parser.add_argument("-z", "--zdegree"  , dest = "zdegree", default = "none", help="Specify z degree")
# parser.add_argument("-a", "--antagonism"  , dest = "antagonism", default = "none", help="Specify antagonism")
# parser.add_argument("-o", "--output"  , dest = "output", default = "none", help="Specify output file name")

# args = parser.parse_args()

# if args.nodes == 'none':
#         Correct_Usage()
# Number_of_Nodes = int(args.nodes)

# if args.trials == 'none':
#         Correct_Usage()
# Number_of_trials = int(args.trials)

# if args.zdegree == 'none':
#         Correct_Usage()
# zdegree = int(args.zdegree)

# if args.antagonism == 'none':
#         Correct_Usage()
# antagonism = float(args.antagonism)

# if args.output == 'none':
#         Correct_Usage()
# output_name = args.output

# print("Nodes               ",Number_of_Nodes)
# print("Number_of_Trials    ",Number_of_trials)
# print("zdegree   ",zdegree)
# print("Output   ",output_name)


#4 possible seed input patterns
def input1(a, b, c): #0,0
    c.vs[a]['state'] = False
    c.vs[b]['state'] = False
    
def input2(a, b, c): #0,1
    c.vs[a]['state'] = False
    c.vs[b]['state'] = True

def input3(a, b, c):#1,0
    c.vs[a]['state'] = True
    c.vs[b]['state'] = False

def input4(a, b, c): #1,1
    c.vs[a]['state'] = True
    c.vs[b]['state'] = True

inputs = [input1, input2, input3, input4]

def generate_network(n_nodes):#not issue
    zvalue = random.randrange(1,10)
    probability = zvalue/(n_nodes-1)
    antagonism = random.uniform(0,1)
    
    G = ig.Graph.Erdos_Renyi(n=n_nodes, p=probability, directed=True)

    ne_antagonistic = int(G.ecount()*antagonism)
    n_antagonistic = int(n_nodes*antagonism)

    G.vs['threshold'] = [random.uniform(0, 1) for _ in range(n_nodes)]
    G.vs['state'] = [False] * n_nodes
    G.es['antagonism'] = [False] * G.ecount()
    G.vs['antagonism'] = [False] * G.vcount()

    
    for ant in random.sample(range(G.ecount()), k=int(ne_antagonistic)):
        G.es[ant]['antagonism'] = True

    for ant in random.sample(range(n_nodes), k=int(n_antagonistic)):
        G.vs[ant]['antagonism'] = True
    
    return G

def calculatefraction(G, n): #not issue
    active = 0
    
    for neighbor in G.neighbors(n, mode="in"):
        edge_id = G.get_eid(neighbor, n)
        if G.vs[neighbor]['state']:
            if G.es[edge_id]['antagonism']:
                active -= 1
            if not G.es[edge_id]['antagonism']:
                active += 1
        
    degree_node = G.degree(n, mode="in")#change mode and loops for directed
    
    #calculate neighbour fraction//
    fraction = active / degree_node if degree_node != 0 else 0

    return fraction

def globalcheck(G): #needs to be done for each input pattern !PROBLEM HERE!
    
    num_labelled = sum(1 for vertex in G.vs if vertex['state'])     
   
    glofraction = num_labelled/(G.vcount())

    if glofraction >= 0.1: #global cascade threshold
        return 1.0
    
    else:
        return 0.0
        

def apply_input_pattern(graph, pattern=inputs):#ok

    nodelist = list(range(graph.vcount()))
    
    random.shuffle(nodelist)
  
    seed1, seed2 = random.sample(nodelist, 2)

    truthtable = np.array([])
    nodetruthtables = np.zeros((len(nodelist), len(pattern)), dtype=int)
    
    for id, input in enumerate(pattern):
        input(seed1, seed2, graph)#ok
        
        for n in nodelist: #individual nodes !issue!
            fraction = calculatefraction(graph, n) #ok

            
            if graph.vs[n]['antagonism'] == False and fraction >= graph.vs[n]['threshold']: #if excitatory
                graph.vs[n]['state'] = True
    
            if graph.vs[n]['antagonism'] and fraction < graph.vs[n]['threshold']:
                graph.vs[n]['state'] = True
                
            else:
                continue

        glo = globalcheck(graph)#ok
        truthtable = np.append(truthtable, glo)
        
        for node in nodelist:
            nodetruthtables[node][id] = int(graph.vs[node]['state'])


        for node in range(graph.vcount()):
            graph.vs[node]['state'] = False
    
    return truthtable, nodetruthtables, graph


def processlogic(trials, nodetables): #evaluates performance  not      
    logiclibrary = {(0.0,0.0,0.0,0.0): 0, (0.0,0.0,0.0,1.0): 1, (0.0,0.0,1.0,0.0): 2, (0.0,0.0,1.0,1.0): 3, 
                    (0.0,1.0,0.0,0.0): 4, (0.0,1.0,0.0,1.0): 5, (0.0,1.0,1.0,0.0): 6, (0.0,1.0,1.0,1.0): 7, 
                    (1.0,0.0,0.0,0.0): 8, (1.0,0.0,0.0,1.0):9, (1.0,0.0,1.0,0.0):10, (1.0,0.0,1.0,1.0):11, 
                    (1.0,1.0,0.0,0.0): 12, (1.0,1.0,0.0,1.0): 13, (1.0,1.0,1.0,0.0): 14, (1.0,1.0,1.0,1.0): 15}

    datalogiclist = []
    nodelogiclist = []
    LTMgfunctions = {str(i): 0 for i in range(16)}
    LTMnfunctions = {str(i): 0 for i in range(16)}

    for array in trials:
        binary = logiclibrary[tuple(array)]
        datalogiclist.append(binary)

    for row in nodetables:
        bin = logiclibrary[tuple(row)]
        nodelogiclist.append(bin)

    # Count occurrences of each logic
    for logic in datalogiclist:
        LTMgfunctions[str(logic)] += 1
    
    for nodelogic in nodelogiclist:
        LTMnfunctions[str(nodelogic)] += 1

    # Create the LTM array
    LTM_array = np.array([LTMgfunctions[str(i)] for i in range(16)])
    node_array = np.array([LTMnfunctions[str(i)] for i in range(16)])

    return LTM_array, node_array 


def namecheckfunc(filename, array):
    while os.path.exists(f"{filename}.npy"):
        digitcheck = any(char.isdigit() for char in filename)

        if not digitcheck:
            filename = f"{filename}_1"

        if digitcheck:
            for char in filename:
                if char.isdigit():
                    digit = int(char)
                    digit += 1
                    filename = f"{filename}_{digit}"
    
    np.save(filename, array)

def namecheckmotif(filename, array):
    while os.path.exists(f"{filename}_motif.npy"):
        digitcheck = any(char.isdigit() for char in filename)

        if not digitcheck:
            filename = f"{filename}_motif1"

        if digitcheck:
            for char in filename:
                if char.isdigit():
                    digit = int(char)
                    digit += 1
                    filename = f"{filename}_motif{digit}"
    
    np.save(f'{filename}_motif', array)



"""
Evolution model - genetic algorithm from scratch

"""
def initialise(n_networks, n_nodes): #generations
    populationnodefunctions = [] #truthtables for each network 
    graphs = []
    count = 0
    while count < n_networks:
        generate = generate_network(n_nodes)
        glo, nodetruth, graph = apply_input_pattern(generate) #truthtable, nodetruthtables,graph

        populationnodefunctions.append(nodetruth)
        graphs.append(graph)
        count += 1

     
    return populationnodefunctions, graphs

def fitness(target, population, graphs): #evaluates performance        
    network_fitness = np.array([])
    for g, nodetable in enumerate(population):
        networkscore = 0
        for node in nodetable:
            nodescore = 0
            for n, output in enumerate(node):
                
                if output == target[n]:
                    nodescore +=1
                if output != target[n]:
                    nodescore -= 0.25
            
            nodefitness = nodescore/len(target)
            networkscore += nodefitness
        
        edgelimit = graphs[g].vcount() * 2.5
        edgecount = graphs[g].ecount()
        #connection penalty
        if edgecount > edgelimit:
            difference = edgecount - edgelimit
            fraction = difference/edgelimit
            networkscore -= fraction

        if networkscore < 0:
            networkscore = 0

        network_fitness = np.append(network_fitness, networkscore)

    return network_fitness

def selection(fitness_score, graphs): #roulette wheel selection - higher fitness has higher share of wheel, chosen at randoms
    norm = np.array([i/np.sum(fitness_score) for i in fitness_score])
    pop_indices = list(range(len(graphs)))
    selectindices = choice(pop_indices, int(len(graphs)/2), p=norm, replace=False)
    selectnetworks = [graphs[index] for index in selectindices]

    return selectnetworks
    
def crossover(selected_population, rate):#needs 0.5 probability, crossover in adjacency matrix or make subgraphs
    
    random.shuffle(selected_population) #select random pair from selected population
    parents = [selected_population[a:a + 2] for a in range(0, len(selected_population),2)]

    for couple in parents:
        if random.random() >= rate and len(couple) == 2:
            crossover_node = random.randint(0, couple[0].vcount())#choose random node
            parent1edges = list((node, neighbor) for node in range(crossover_node, couple[0].vcount()) for neighbor in couple[0].neighbors(node)) #make list of edges from node crossover point
            t1 = couple[0].vs['threshold']#need to get both edge and node attributes
            n1 = couple[0].vs['antagonism']
            e1 = couple[0].es['antagonism']
            thresholds1 = [t1[i] for i in range(crossover_node, couple[0].vcount())]
            nodetypes1 = [n1[i] for i in range(crossover_node, couple[0].vcount())]
            edgetypes1 = []
            
            for node, neighbor in parent1edges:
            # Get the edge ID for the current edge
                if couple[0].are_adjacent(node, neighbor):
                    edge_id = couple[0].get_eid(node, neighbor)
                    edgetypes1.append(e1[edge_id])

            
            
            parent2edges = list((node, neighbor) for node in range(crossover_node, couple[1].vcount()) for neighbor in couple[1].neighbors(node))
            t2 = couple[1].vs['threshold']#need to get both edge and node attributes
            n2 = couple[1].vs['antagonism']
            e2 = couple[1].es['antagonism']

            thresholds2 = [t2[i] for i in range(crossover_node, couple[1].vcount())]
            nodetypes2 = [n2[i] for i in range(crossover_node, couple[1].vcount())]
            edgetypes2 = []
            
            for node, neighbor in parent2edges:
            # Get the edge ID for the current edge
                if couple[1].are_adjacent(node, neighbor):
                    edge_id = couple[1].get_eid(node, neighbor)
                    edgetypes2.append(e2[edge_id])

           
            #parent1 crossover
            edge_ids_to_remove1 = [couple[0].get_eid(edge[0], edge[1]) for edge in parent1edges if couple[0].are_adjacent(edge[0], edge[1])]
            couple[0].delete_edges(edge_ids_to_remove1)#remove original edges
            
            # edge_ids_to_add1 = [couple[1].get_eid(edge[0], edge[1]) for edge in parent2edges if couple[1].are_connected(edge[0], edge[1])]
            couple[0].add_edges(parent2edges)#replace w new edges
            
            for i, y in enumerate(range(crossover_node, couple[1].vcount())):
                couple[0].vs[y]['threshold'] = thresholds2[i]
                couple[0].vs[y]['antagonism'] = nodetypes2[i]
            
            for j, x in enumerate(parent2edges):
                edgeid = couple[0].get_eid(x[0], x[1])
                try:
                    couple[0].es[edgeid]['antagonism'] = edgetypes2[j]
                
                except IndexError:
                    couple[0].es[edgeid]['antagonism'] = False

                else:
                    couple[0].es[edgeid]['antagonism'] = edgetypes2[j]
            
            
            #parent2 crossover
            edge_ids_to_remove2 = [couple[1].get_eid(edge[0], edge[1]) for edge in parent2edges if couple[1].are_adjacent(edge[0], edge[1])]
            couple[1].delete_edges(edge_ids_to_remove2)#remove original edges
            
            # edge_ids_to_add2 = [couple[0].get_eid(edge[0], edge[1]) for edge in parent1edges if couple[0].are_connected(edge[0], edge[1])]
            couple[1].add_edges(parent1edges)#replace w new edges
            
            for i, y in enumerate(range(crossover_node, couple[0].vcount())):
                couple[1].vs[y]['threshold'] = thresholds1[i]
                couple[1].vs[y]['antagonism'] = nodetypes1[i]
            
            for j, x in enumerate(parent1edges):
                edgeid = couple[1].get_eid(x[0], x[1])
                try:
                    couple[1].es[edgeid]['antagonism'] = edgetypes1[j]
                
                except IndexError:
                    couple[1].es[edgeid]['antagonism'] = False

                else:
                    couple[1].es[edgeid]['antagonism'] = edgetypes1[j]
                
            for edges in range(couple[0].ecount()):
                if couple[0].es[edges]['antagonism'] == None:
                    couple[0].es[edges]['antagonism'] = False

            for edges in range(couple[1].ecount()):
                if couple[1].es[edges]['antagonism'] == None:
                    couple[1].es[edges]['antagonism'] = False

    children = [cross for couple in parents for cross in couple ]
    
    return children

def mutate(children, mut_rate):#need to apply LTM here

    new_gen = [child for child in children for _ in range(2)]
    
    mutated = [] #actual mutant graphs
    mut_population = [] #truth tables of each network output
    types = ['type', 'threshold', 'edge']
    
    for network in new_gen:
        if random.random() >= mut_rate:
            typechosen = random.choice(types)
            node = random.choice(list(range(network.vcount())))
            edge = random.choice(list(range(network.ecount())))
            
            if typechosen == 'type': #choose a random node and then random choice of excitatory or antagonistic
                network.es[edge]['antagonism'] = random.choice([True, False])
            
            if typechosen == 'threshold':
                network.vs[node]['threshold'] = random.uniform(0, 1)
            
            if typechosen == 'edge': #remove or add edge at random node
                neighbours = list(network.neighbors(node))
                if neighbours:
                    choose_edge = random.choice(neighbours)
                    edgechoice = random.choice(['remove', 'add'])
                    if network.are_adjacent(node, choose_edge):
                        edgeid = network.get_eid(node, choose_edge)
                    else:
                        continue

                    if edgechoice == 'remove':
                        if network.are_adjacent(node, choose_edge):
                            network.delete_edges(edgeid)
            
                    if edgechoice == 'add':
                        if network.are_adjacent(node, choose_edge):
                            network.add_edge(node, choose_edge)

        muttruth, mutnode, mutant = apply_input_pattern(network)
        mut_population.append(mutnode)
        mutated.append(mutant)
    

    return muttruth, mut_population, mutated



def GA(n_networks, nodes, target, n_gen): #call this function to run GA
    
    population, graphs = initialise(n_networks,nodes)
    cascadetables = []
    
    for gen in range(n_gen):
        fitness_score = fitness(target, population, graphs)
        selected = selection(fitness_score, graphs)
        children = crossover(selected, 0.5)
        cascade, population, graphs = mutate(children, 0.5)
        cascadetables.append(cascade)
    
    return population, cascadetables




