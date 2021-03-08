import numpy as np
import matplotlib.pyplot as plt
import time
from evolving_mechanism import *
from cdlib import algorithms

import igraph
import entropy as ep
import igraph as ig


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
    
    t0 =time.time()
    #a = g.pagerank(vertices=[0,1])
    #a = g.betweenness(vertices=[3,10])
    a = g.eigenvector_centrality()
    print("time is",time.time()-t0)
    #print(a)
    # print(t1-t0,t2-t1,t3-t2)
    # print(igraph.summary(g))
    # a,b = nx.cn_soundarajan_hopcroft(G.to_undirected(), [(0, 2),(1,2)])
    # c,d = nx.cn_soundarajan_hopcroft(G.to_undirected(), [(0, 2),(1,2)])
    # # a,b = cn_soundarajan_hopcroft(G, ebunch=None, community='community')
    # print(a,b,c,d)
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
    
    uG.nodes[1]["community"] =1
    uG_c = uG.copy()
    uG.nodes[1]["community"] = 2
    print(uG_c.nodes[1])
        
    

        

    # print(igraph.summary(g))
    # print(g)
    
    # print(g.degree(2))
    
    
    # print(a[2],b[2])
    # print(g.nodes,type(g.nodes()))

def test_index():
    el = EdgeList()
    G = el.read_BA(1.0,10)
    
    # a,b= nx.resource_allocation_index(g.to_undirected(),[(0,1),(1,2)])#done
    # a,b= nx.jaccard_coefficient(g.to_undirected(),[(0,1),(1,2)]) #OK
    # a,b=nx.adamic_adar_index(g.to_undirected(),[(0,1),(1,2)]) #ok
    # preferential_attachment
    # common_neighbor_centrality

    # print(g.vs)
    # g.add_edges(list(G.edges()))
    

    
    
    nx.info(G)
    uG = G.to_undirected()
    uG.add_edge(104,105)
    
    g = igraph.Graph()
    g.add_vertices(list(uG.nodes()))
    for e in uG.edges():
        # print(ver,ver.index,ver["name"],type(ver["name"]))
        s = g.vs.select(name = e[0])[0].index
        d = g.vs.select(name = e[1])[0].index
        g.add_edges([(s,d)])
        pass 
    
    #eigenvector
    print(g.evcent())
    print(list(nx.eigenvector_centrality(uG).values()))
     
    #communicability_exp
    #c = nx.communicability_exp(uG)
    # print(c[0][1], c[1][0])
    
    #laplacian
    # lap = g.laplacian()
    # print(True in np.isnan(np.array(lap)))
    
    
    # eigenvector
    #hits
    # print(nx.hits(uG)[0])
    # print(nx.hits(uG)[1])
    
    
    # s = g.vs.select(name = 9)[0].index
    # d = g.vs.select(name = 104)[0].index
    # L = np.linalg.pinv((g.laplacian()))
    # print(L[s,d])
    
    # print(g.personalized_pagerank(reset_vertices=2))
    # print(g.personalized_pagerank(reset_vertices=1))
    # print(g.pagerank())
    # print(nx.pagerank(uG,max_iter=1,tol=0.1))

    # n = len(g.laplacian())
    # mfi = np.linalg.inv(np.identity(n)+np.array(g.laplacian()))
    # print(mfi[9,10],mfi[1,0])
    
    
    

    # print(uG.nodes)
    # L = nx.laplacian_matrix(uG)
    # print(L)
    # L_pinv = np.linalg.pinv(L.toarray())
    # print(L_pinv)
    # print(list(uG.nodes()).index(104))
    
    # loca = index_local()
    # a, b =loca.average_commute_time(uG, 8,9,104)
    # print(a,b)
    
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
    # a = index_global._h_index([])
    
    # print(a)
    
   
    
    # a=nx.adamic_adar_index(G,[(1,1),(1,1)]) #ok
    # print(G.nodes)
    # a  = list(G.nodes)
    # a.remove(0)
    # print(a)
    # print(G.nodes)
    # a,b = nx.common_neighbor_centrality(G.to_undirected(),[(0,1),(1,2)])
    # print(a[2],b[2])
    
    
    
    pass


test_index()
# test_networkx()



#test_index()

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



from multiprocessing.pool import ThreadPool as Pool
import os, time, random

d = 2
lis = []
def long_time_task(name):
    print('Run task %s (%s)...' % (name, os.getpid()))
    start = time.time()
    t = random.random() * 3
    time.sleep(t)
    end = time.time()
    global lis
    lis.append(t)
    print('Task %s runs %0.2f seconds.' % (name, (end - start)))
    return 1,d, d**2
    


    
def test_time():
    d = {'indegree': 0.024885177612304688, 'degree': 0.009556055068969727, 'continuous_degree': 0.25249171257019043, 'k1sum': 0.030849456787109375, 'k2sum': 0.04731035232543945, 'closeness_uG': 0.6424429416656494, 'harmonic_uG': 0.6728346347808838, 'pagerank_ig': 0.36860203742980957, 'h_index_uG': 0.10823535919189453, 'kshell_uG': 1.062082290649414, 'resource_allocation_uG': 0.10199880599975586, 'jaccard_coefficient_uG': 0.07625079154968262, 'adamic_adar_index_uG': 0.0676875114440918, 'preferential_attachment_uG': 0.020041227340698242, 'cn_soundarajan_hopcroft_uG': 0.07007884979248047, 'ra_index_soundarajan_hopcroft_uG': 0.025110483169555664, 'within_inter_cluster_uG': 0.02451324462890625, 'common_neighbor_centrality_uG': 0.09677505493164062, 'efficiency_uG': 0.03606915473937988, 'distance_uG': 0.0272519588470459, 'distance_cutoff_uG': 0.02762126922607422, 'community': 0.9620838165283203}
    s = np.sum(list(d.values()))
    for i in d:
        d[i] = d[i]/s
    # print(d)
    
    a = sorted(d.items(), key=lambda x: x[1], reverse=True)
    print(a)  
    import multiprocessing
    cores = multiprocessing.cpu_count()
    print(cores)


def output_text():
    filename = "./experimental_results/output.txt"
    with open(filename, 'a') as f:
        print(var, file=f)
    pass

def np_test():
    pass
    # el = EdgeList()
    # file_name = "./experimental_results/hepph_all.xlsx"
    # el.load_records(file_name)
    # el.cut_and_smooth_normalize_records(normalize_indices=["degree"],smooth_length=300)
    # print(el.records["degree"])
    # print(el.records["degree"][0:1000])
    
    # m = 6259
    # slice_num = int(m/10)
    # for d in range(10):
    #     if d != 9:
    #         print(d,d*slice_num,(d+1)*slice_num)
    #     else:
    #         print(d,d*slice_num,m-1)
        
def test_slice():
    nums = [1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5]
    print(ep.slice(nums,4))
    #print(list(range(0,9)))
    

    
    
