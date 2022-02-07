import collections
import matplotlib.pyplot as plt
import random
import utils
import sys
sequence_reads, qualities = utils.read_fastq('TeleTubby.fastq')
print("number of reads: ", len(sequence_reads))

def build_graph(k_mers):
    edges = []
    nodes = set()
    cycleDict = {} #dictionary to keep track of prefix/suffix
    k_merLength = len(k_mers[0]) #storing length to reduce operations
    
    for k_mer in k_mers: 
        pre = k_mer[:k_merLength-1] #create prefix key
        suf = k_mer[1:] #create suffix
        if pre in cycleDict.keys():
            cycleDict.get(pre).add(suf) #if the prefix already exists, add another word/suffix it connects to
        else:
            cycleDict[pre]= set() #if the prefix doesn't already exist, add the word/suffix
            cycleDict[pre].add(suf)
    
    #create list of nodes and edges for graph viz
    for pre in cycleDict.keys(): #for every prefix
        nodes.add(pre)
        for suf in cycleDict.get(pre):
            nodes.add(suf)
            edges.append((pre,suf))
    #calculate indegree and outdegree of each node
    degree = {i:[0,0] for i in nodes}
    for pre in cycleDict.keys():
        for suf in cycleDict[pre]:
            degree[suf][0] += 1
            degree[pre][1] += 1
            
    for pre in cycleDict.keys():
        cycleDict[pre] = list(cycleDict[pre])
    #find the unbalanced nodes
    unbalanced = [(i,degree[i][1]-degree[i][0]) for i in degree.keys() if degree[i][1]-degree[i][0]]
    unbalanced.sort(key = lambda x: -x[1])
    print('unbalanced nodes: ', unbalanced)
    print("number of nodes: ", len(nodes))
    return nodes, edges, cycleDict, unbalanced


def fast_euler(graph, cur):
    path = []
    stack = []
    stack.append(cur)
    while stack:
        #print("stack: ", stack)
        #print(graph)
        if graph[cur]:
            stack.append(cur)
            adj = graph[cur][-1]
            # print(graph[cur])
            graph[cur].remove(adj)
            cur =adj
        else:
            path.append(cur)
            cur = stack.pop()
    return list(reversed(path))
# graph = {'a':['b'], 'b':['c'], 'c':['a']}
# print(fast_euler(graph,'a'))

nodes, edges, cycleDict, unbalanced = build_graph(sequence_reads)
start = unbalanced[0][0]
cycleDict[unbalanced[1][0]] = [start]
path = fast_euler(cycleDict,start)
sequence = path[0]+ "".join([i[-1] for i in path[1:-1]])
print(sequence, len(sequence))

sequence_reads_covid, qualities_covid = utils.read_fastq('SARS-CoV2.fastq')
nodes_covid, edges_covid, cycle_dict, unbalanced = build_graph(sequence_reads_covid)

# Your code here
fd = {}
for read in sequence_reads_covid:
    if read in fd:
        fd[read]+=1
    else:
        fd[read]=1

top = sorted(list(fd.keys()), key=lambda x: -fd[x])
print([(i,fd[i]) for i in top[:3]])

start = unbalanced[0][0]
cycle_dict[unbalanced[1][0]] = [start]
path = fast_euler(cycle_dict,start)
sequence = path[0]+ "".join([i[-1] for i in path[1:-1]])
print(sequence, len(sequence))