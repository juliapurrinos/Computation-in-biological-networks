#!/usr/bin/env python
import numpy as np
import igraph as ig
import random
import sys, argparse
import os


'''
activation iteration model: 
iterates through each possible network and for each possible network it iterates through 
each possible input pattern then saves the data in a numpy array.
Each for loop of network generation = 1 trial
'''

def Correct_Usage():
        print("Correct Way to run program is:")
        print("progname -n number of nodes -t number of trials -z zvalue -a antagonism -o outputname")
        sys.exit(1)


# ==== Main program starts here

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nodes"    , dest = "nodes",   default = "none", help="Nodes")
parser.add_argument("-t", "--trials"  , dest = "trials", default = "none", help="Trials")
parser.add_argument("-z", "--zdegree"  , dest = "zdegree", default = "none", help="Specify z degree")
parser.add_argument("-a", "--antagonism"  , dest = "antagonism", default = "none", help="Specify antagonism")
parser.add_argument("-o", "--output"  , dest = "output", default = "none", help="Specify output file name")

args = parser.parse_args()

if args.nodes == 'none':
        Correct_Usage()
Number_of_Nodes = int(args.nodes)

if args.trials == 'none':
        Correct_Usage()
Number_of_trials = int(args.trials)

if args.zdegree == 'none':
        Correct_Usage()
zdegree = int(args.zdegree)

if args.antagonism == 'none':
        Correct_Usage()
antagonism = float(args.antagonism)

if args.output == 'none':
        Correct_Usage()
output_name = args.output

print("Nodes               ",Number_of_Nodes)
print("Number_of_Trials    ",Number_of_trials)
print("zdegree   ",zdegree)
print("Output   ",output_name)


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

def generate_network(n_nodes, zvalue, antagonism):#not issue
    probability = zvalue/(n_nodes-1)
    
    G = ig.Graph.Erdos_Renyi(n=n_nodes, p=probability, directed=True)

    ne_antagonistic = int(G.ecount()*antagonism)
    n_antagonistic = int(n_nodes*antagonism)

    G.vs['threshold'] = [random.uniform(0, 1) for _ in range(n_nodes)]
    G.vs['state'] = [False] * n_nodes
    G.vs['antagonism'] = [False] * G.vcount()
    G.es['antagonism'] = [False] * G.ecount()

    
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
        

def apply_input_pattern(graph, pattern):#ok

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
    
    return truthtable, nodetruthtables


def processlogic(glob, nodetable): #evaluates performance  not      
    logiclibrary = {(0.0,0.0,0.0,0.0): 0, (0.0,0.0,0.0,1.0): 1, (0.0,0.0,1.0,0.0): 2, (0.0,0.0,1.0,1.0): 3, 
                    (0.0,1.0,0.0,0.0): 4, (0.0,1.0,0.0,1.0): 5, (0.0,1.0,1.0,0.0): 6, (0.0,1.0,1.0,1.0): 7, 
                    (1.0,0.0,0.0,0.0): 8, (1.0,0.0,0.0,1.0):9, (1.0,0.0,1.0,0.0):10, (1.0,0.0,1.0,1.0):11, 
                    (1.0,1.0,0.0,0.0): 12, (1.0,1.0,0.0,1.0): 13, (1.0,1.0,1.0,0.0): 14, (1.0,1.0,1.0,1.0): 15}

    LTMgfunctions = {str(i): 0 for i in range(16)}
    LTMnfunctions = {str(i): 0 for i in range(16)}

    for row in nodetable:
        nodelogic = logiclibrary[tuple(row)]
        LTMnfunctions[str(nodelogic)] += 1
    
    glologic = logiclibrary[tuple(glob)]
    LTMgfunctions[str(glologic)] += 1

    global_array = np.array([LTMgfunctions[str(i)] for i in range(16)])
    node_array = np.array([LTMnfunctions[str(i)] for i in range(16)])
    
    return global_array, node_array

def namecheckglo(filename, array):
    while os.path.exists(f"{filename}_glo.npy"):
        digitcheck = any(char.isdigit() for char in filename)

        if not digitcheck:
            filename = f"{filename}_1"

        if digitcheck:
            for char in filename:
                if char.isdigit():
                    digit = int(char)
                    digit += 1
                    filename = f"{filename}_{digit}"
    
    np.save(f'{filename}_glo', array)

def namechecknode(filename, array):
    while os.path.exists(f"{filename}_node.npy"):
        digitcheck = any(char.isdigit() for char in filename)

        if not digitcheck:
            filename = f"{filename}_1"

        if digitcheck:
            for char in filename:
                if char.isdigit():
                    digit = int(char)
                    digit += 1
                    filename = f"{filename}_{digit}"
    
    np.save(f'{filename}_node', array)

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


def run_trial(n_trials, n_nodes, zvalue, antagonism, filename):
    trials = np.empty((0, 16), dtype=int)
    nodestack = np.empty((0, 16), dtype=int)
    motifcount = np.empty((0, 16), dtype=int)
    
    for _ in range(n_trials):
        generate = generate_network(n_nodes, zvalue, antagonism)#ok
        graphoutput, nodeoutput = apply_input_pattern(generate, inputs)#ok
        glocount, nodecount = processlogic(graphoutput, nodeoutput)

        trials = np.vstack((trials, glocount))
        nodestack = np.vstack((nodestack, nodecount))

        motifs = np.array([generate.motifs_randesu(size=3)])
        motifcount = np.vstack((motifcount, motifs))

    namechecknode(filename, nodestack)
    namecheckglo(filename, trials)
    namecheckmotif(filename, motifcount)

run_trial(Number_of_trials, Number_of_Nodes, zdegree, antagonism, output_name)