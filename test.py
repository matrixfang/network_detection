import numpy as np
import matplotlib.pyplot as plt
import time
from evolving_mechanism import *
from cdlib import algorithms

import igraph



def test_container():
    pass
    x = list(np.random.rand(1000) -1.0/500*np.arange(0,1000))
    y = np.ones(1000)
    #print(np.arange(0,1000))
    N=200
    weights = np.ones(N)/N

    sma = np.convolve(x,weights,'same')

    t = np.arange(0,len(x))
    print(len(sma),len(x),len(t))
    plt.plot(t[N-1:-N+1],x[N-1:-N+1])
    plt.plot(t[N-1:-N+1],sma[N-1:-N+1])
    print(len(t[N-1:-N+1]))
    print(t[N-1:])
    plt.show()


    a = np.array([5,6,7,8,1])

    a.sort()
    from collections import deque

    d = deque()
    d.append(1)
    d.append('a')
    d.append(1)
    print(d)

    t0 = time.time()
    i = 0
    dd = []
    for j in range(100000):
        dd.append((i,i+1,i+2))
        
        #print(i)
    print("time is ", time.time()-t0)


    t0 = time.time()
    i = 0
    import array
    arr0 = array.array('d', [1])
    arr1 = array.array('d', [2])
    arr2 = array.array('d', [3])
    for j in range(100000):
        arr0.append(i)
        arr1.append(i+1)
        arr2.append(i+2)
        
        #print(i)
    print("time is ", time.time()-t0)
    
    
    
def test_scale_free():
    el = EdgeList()
    g = el.read_BA(1.0)
    degree_sequence = sorted([d for n,d in g.degree()],reverse=True)
    import powerlaw
    fit = powerlaw.Fit(degree_sequence)
    #print()
    fig2 = fit.plot_pdf(color='b',linewidth=2)
    fit.power_law.plot_pdf(color ='g',linestyle='--',ax=fig2)
    plt.annotate("$\gamma=$"+str(fit.power_law.alpha),xy=(0.6,0.6),
                 xytext=(0.6,0.6),xycoords='axes fraction')
    # pos = nx.spring_layout(g)
    
    # nx.draw_networkx_nodes(g, pos, node_size=20)
    # nx.draw_networkx_edges(g, pos, alpha=0.4)
    plt.show()
    
    
def python3_test():
    d = {"a":1,"b":2}
    print(d.values(),type(d.values()))
    print(d.keys(),type(d.keys()))
    print(list(d.keys()))
    for i in d.items():
        print(type(i))


def add_attribute_community():
    # G = nx.karate_club_graph()
    el = EdgeList()
    
    t0 = time.time()
    G = el.read_BA(1.0,10000)
    print(nx.info(G))
    print("end of networkx information ")
    
    g = igraph.Graph(directed=False)

    
    t1 = time.time()
    
    # print(G.nodes())
    # print(G.edges())
    g.add_vertices(list(G.nodes()))
    g.add_edges(list(G.edges()))
    # g = igraph.Graph.from_networkx(G)
    p_louvain = g.community_multilevel()
    # print(G.nodes())
    # print(p_louvain)
    mod = g.modularity(p_louvain)
    print("modularity is ", mod)
    
    t2 = time.time()
    G.to_undirected()
    t3 = time.time()
    cluster_index = 0
    for cluster in p_louvain:
        cluster_index += 1
        for node in cluster:
            G.nodes[node]["community"] = cluster_index
    

    print(t1-t0,t2-t1,t3-t2)
    print(igraph.summary(g))
    a,b = nx.cn_soundarajan_hopcroft(G.to_undirected(), [(0, 2),(1,2)])
    c,d = nx.cn_soundarajan_hopcroft(G.to_undirected(), [(0, 2),(1,2)])
    # a,b = cn_soundarajan_hopcroft(G, ebunch=None, community='community')
    print(a,b,c,d)
    # print(coms.communities)
    pass


def test_networkx():
    el = EdgeList()
    G = el.read_BA(1.0,1000)
    # a,b= nx.resource_allocation_index(g.to_undirected(),[(0,1),(1,2)])#done
    # a,b= nx.jaccard_coefficient(g.to_undirected(),[(0,1),(1,2)]) #OK
    # a,b=nx.adamic_adar_index(g.to_undirected(),[(0,1),(1,2)]) #ok
    # preferential_attachment
    # common_neighbor_centrality
    # file_name = './real_evolving_networks/sx-mathoverflow-a2q.txt'
    # el = EdgeList() #build edgelist object
    # index4normalize = ["indegree","degree","katz","k1sum","k2sum","resource_allocation"]
    # index4cutoff = ['distance','distance_cutoff']
    # el.load_index_for_recorded(index4normalize+index4cutoff)
    
    # el.read_edgelist_file(file_name)
    # G = el.evolve(max_edge_number = 10000)

    # nx.info(G)
    
    t0 = time.time()
    uG = G.to_undirected()
    t1 =time.time()
    print("ug time is",t1-t0)

    
    g = igraph.Graph()
    g.add_vertices(list(G.nodes()))
    # print(g.vs)
    # g.add_edges(list(G.edges()))
    
    for e in G.edges():
        # print(ver,ver.index,ver["name"],type(ver["name"]))
        s = g.vs.select(name = e[0])[0].index
        d = g.vs.select(name = e[1])[0].index
        g.add_edges([(s,d)])
        pass
        
    

        

    # print(igraph.summary(g))
    # print(g)
    
    # print(g.degree(2))
    
    
    # print(a[2],b[2])
    # print(g.nodes,type(g.nodes()))

def test_index():
    el = EdgeList()
    G = el.read_BA(1.0,100)
    # a,b= nx.resource_allocation_index(g.to_undirected(),[(0,1),(1,2)])#done
    # a,b= nx.jaccard_coefficient(g.to_undirected(),[(0,1),(1,2)]) #OK
    # a,b=nx.adamic_adar_index(g.to_undirected(),[(0,1),(1,2)]) #ok
    # preferential_attachment
    # common_neighbor_centrality
    
    # g = nx.Graph()
    # g.add_edge(1,1)
    
    nx.info(G)
    uG = G.to_undirected()
    # a = nx.closeness_centrality(uG,u=0)
    # d = nx.betweenness_centrality(uG)
    # d = nx.current_flow_betweenness_centrality(uG)
    # a = nx.harmonic_centrality(uG,[1,2])
    # a =  nx.rich_club_coefficient(uG)
    # a = nx.pagerank_numpy(uG)
    # print(a,type(a))
    
    # d = nx.katz_centrality_numpy(G)
    # d = nx.eigenvector_centrality_numpy(uG)
    # d = nx.closeness_centrality(uG)
    #uG.add_edge(102,103)
   # d = nx.current_flow_betweenness_centrality(uG)
    a = index_global._h_index([])
    
    print(a)
    
   
    
    # a=nx.adamic_adar_index(G,[(1,1),(1,1)]) #ok
    # print(G.nodes)
    # a  = list(G.nodes)
    # a.remove(0)
    # print(a)
    # print(G.nodes)
    # a,b = nx.common_neighbor_centrality(G.to_undirected(),[(0,1),(1,2)])
    # print(a[2],b[2])
    
    
    
    pass

# test_index()

# test_networkx()


# add_attribute_community()

test_index()

def function_test():
    s = "dfasdfasdfa"
    print(s[-2:])

# function_test()

def h_index():
    from collections import Counter
    indexList = [1,1,1,2,2,2,2,5,5,5,5,5,6]
    HCounter = Counter(indexList)
    CounterKeys = [i for i in HCounter]
    CounterKeys = sorted(CounterKeys,reverse=True) 
    
    CounterValues = [HCounter[i] for i in CounterKeys]
    print(CounterKeys,CounterValues)
    for index in range(0,len(CounterValues)):
        if CounterKeys[index]<=sum(CounterValues[0:index+1]):
            break
    print(CounterKeys[index])
    return CounterKeys[index]
    
    pass


