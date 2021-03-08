import networkx as nx
import numpy as np
import random 
import matplotlib.pyplot as plt
import entropy as ep
import pandas as pd
from npeet import entropy_estimators as ee
import os
from tqdm import tqdm
import igraph
from collections import Counter
# from multiprocessing import Pool
from multiprocessing.pool import ThreadPool as Pool
import multiprocessing as mp
import time
import scipy as sp

def sent_notifer():
    command = "terminal-notifier -activate 'com.microsoft.VSCode' -sender 'com.microsoft.VSCode' -title 'python notifer' -message 'Your python program has finished!' "
    os.system(command)
    os.system('say "your program has finished"')
    # com.apple.Safari
    
    pass

    


class index_local(object):
    def __init__(self):
        self.laplacian_matrix_pinv = None
        pass
    def pesodoinv_laplacian(self):
        pass
    @staticmethod
    def test():
        pass
    @staticmethod
    def _cutoff_distance(g,ni,nj):
        if nj in g.neighbors(ni):
            return 1
        for neigh in g.neighbors(ni):
            if nj in g.neighbors(neigh):
                return 2
            else:
                pass
        return 3
   
    @staticmethod 
    def cutoff_distance(uG,ni,nj,rand_node):
        return index_local._cutoff_distance(uG,ni,nj),index_local._cutoff_distance(uG,ni,rand_node)
        
    @staticmethod
    def _shortest_path_length(uG,ni,nj):
        try:
            length = nx.shortest_path_length(uG, source=ni, target=nj)
        except nx.exception.NetworkXNoPath:
            return 100
        else:
            return length
    @staticmethod
    def shortest_path_length(uG,ni,nj,rand_node):
        return index_local._shortest_path_length(uG,ni,nj),index_local._shortest_path_length(uG,ni,rand_node)
    
    @staticmethod
    def resource_allocation(uG,ni,nj,rand_node):
        a,b = nx.resource_allocation_index(uG,[(ni,nj),(ni,rand_node)])
        return a[2],b[2]
    @staticmethod 
    def jaccard_coefficient(uG, ni,nj,rand_node):
        a,b = nx.jaccard_coefficient(uG,[(ni,nj),(ni,rand_node)])
        return a[2],b[2]
    @staticmethod
    def adamic_adar_index(uG,ni,nj,rand_node):
        a,b = nx.adamic_adar_index(uG, [(ni, nj), (ni, rand_node)])
        return a[2],b[2]
    @staticmethod
    def preferential_attachment(uG, ni,nj,rand_node):
        a,b = nx.preferential_attachment(uG,[(ni,nj),(ni,rand_node)])
        return a[2],b[2]
    
    @staticmethod
    def cn_soundarajan_hopcroft(uG,ni,nj,rand_node):
        a,b = nx.cn_soundarajan_hopcroft(uG,[(ni,nj),(ni,rand_node)])
        return a[2],b[2]
    @staticmethod
    def ra_index_soundarajan_hopcroft(uG, ni, nj,rand_node):
        a,b = nx.ra_index_soundarajan_hopcroft(uG,[(ni,nj),(ni,rand_node)])
        return a[2],b[2]
    @staticmethod
    def within_inter_cluster(uG,ni,nj,rand_node):
        a,b = nx.within_inter_cluster(uG,[(ni,nj),(ni,rand_node)])
        return a[2],b[2]
    @staticmethod
    def common_neighbor_centrality(uG,ni,nj,rand_node):
        alpha = 0.8
        a =  0.8*len(list(nx.common_neighbors(uG,ni,nj)))+(1-alpha)*len(uG.nodes())/index_local._shortest_path_length(uG,ni,nj)
        b =  0.8*len(list(nx.common_neighbors(uG,ni,rand_node))) + (1-alpha)*len(uG.nodes())/index_local._shortest_path_length(uG,ni,rand_node)
        
        # a,b = nx.common_neighbor_centrality(uG,[(ni,nj),(ni,rand_node)])
        return a,b
    
    @staticmethod
    def salton(uG, ni, nj, rand_node):
        a = len(list(nx.common_neighbors(uG,ni,nj)))/(np.sqrt(uG.degree(nj) * uG.degree(ni)))
        b = len(list(nx.common_neighbors(uG,ni,nj)))/(np.sqrt(uG.degree(ni) * uG.degree(rand_node)))
        return a,b
    @staticmethod
    def sorensen(uG, ni, nj, rand_node):
        a = 2*len(list(nx.common_neighbors(uG,ni,nj)))/(uG.degree(nj)+uG.degree(ni))
        b = 2*len(list(nx.common_neighbors(uG,ni,nj)))/(uG.degree(ni)+uG.degree(rand_node))
        return a, b
    @staticmethod
    def HPI(uG, ni, nj, rand_node):
        a = len(list(nx.common_neighbors(uG,ni,nj)))/min(uG.degree(ni),uG.degree(nj))
        b = len(list(nx.common_neighbors(uG,ni,nj)))/min(uG.degree(ni),uG.degree(rand_node))  
        return a, b
    @staticmethod
    def HDI(uG, ni, nj, rand_node):
        a = len(list(nx.common_neighbors(uG,ni,nj)))/max(uG.degree(ni),uG.degree(nj))
        b = len(list(nx.common_neighbors(uG,ni,nj)))/max(uG.degree(ni),uG.degree(rand_node))  
        return a, b
    @staticmethod
    def LHN1(uG, ni, nj, rand_node):
        a = len(list(nx.common_neighbors(uG,ni,nj)))/(uG.degree(ni)*uG.degree(nj))
        b = len(list(nx.common_neighbors(uG,ni,nj)))/(uG.degree(ni)*uG.degree(rand_node))  
        return a, b

    def average_commute_time(self,g, ni,nj, rand_node):
        lap = np.array(g.laplacian())
        try:
            l = sp.linalg.pinv(lap)
            
        except np.linalg.LinAlgError:
            lap = g.laplacian()
            print("matrix is lap")
            print(True in np.isnan(lap))
        
        self.laplacian_matrix_pinv  = l 
        ni_index = g.vs.select(name = ni)[0].index
        nj_index = g.vs.select(name = nj)[0].index
        rand_index = g.vs.select(name = rand_node)[0].index

        a = l[ni_index,ni_index]+l[nj_index,nj_index]-2*l[ni_index,nj_index]
        b = l[ni_index,rand_index]+l[nj_index,rand_index]-2*l[ni_index,rand_index]
        return a,b

    def cos(self,g,ni,nj,rand_node):
        l = self.laplacian_matrix_pinv 
        ni_index = g.vs.select(name = ni)[0].index
        nj_index = g.vs.select(name = nj)[0].index
        rand_index = g.vs.select(name = rand_node)[0].index
        
        a = l[ni_index,nj_index]/np.sqrt(abs(l[ni_index,ni_index]*l[nj_index,nj_index]))
        b = l[ni_index,rand_index]/np.sqrt(abs(l[ni_index,ni_index]*l[rand_index,rand_index]))
        self.laplacian_matrix_pinv  = None
        return a,b
    
    @staticmethod
    def random_walk_with_restart(g,ni,nj,rand_node):
        ni_index = g.vs.select(name = ni)[0].index
        nj_index = g.vs.select(name = nj)[0].index
        rand_index = g.vs.select(name = rand_node)[0].index
        a = g.personalized_pagerank(reset_vertices=ni_index)[nj_index] + g.personalized_pagerank(reset_vertices=nj_index)[ni_index]
        b = g.personalized_pagerank(reset_vertices=ni_index)[rand_index] + g.personalized_pagerank(reset_vertices=rand_index)[ni_index] 
        
        return a, b
    @staticmethod
    def matrix_forest_index(g, ni,nj,rand_node):
        n = len(g.vs)
        mfi = np.linalg.pinv(np.identity(n)+np.array(g.laplacian()))
        ni_index = g.vs.select(name = ni)[0].index
        nj_index = g.vs.select(name = nj)[0].index
        rand_index = g.vs.select(name = rand_node)[0].index
        return mfi[ni_index,nj_index],mfi[ni_index,rand_index]
        
    @staticmethod
    def efficiency(uG,ni,nj,rand_node):
        a = nx.efficiency(uG,ni,nj)
        b = nx.efficiency(uG,ni,rand_node)
        return a,b
    
    # bascially all the link prediction index is locally index, can be viewed as some kind of similarity
    

class index_global(object):
    def __init__(self):
        # ni does not affect the value of index_global
        pass
    @staticmethod
    def _h_index(indexList):
        """
        output h-index of the indexLis looks like [6,5,5,5,5,2,2,2,2,2,1,1,1,1,1,1]
        """
        HCounter = Counter(indexList)
        CounterKeys = [i for i in HCounter]
        CounterKeys = sorted(CounterKeys,reverse=True) 
        
        CounterValues = [HCounter[i] for i in CounterKeys]
        #print(CounterKeys,CounterValues)
        for index in range(0,len(CounterValues)):
            if CounterKeys[index]<=sum(CounterValues[0:index+1]):
                break
        return CounterKeys[index]
    @staticmethod
    def _ink1sum(g,ni,nj):
        num = 0
        for neigh in g.predecessors(nj):
            num += g.in_degree(neigh)
        return num
    @staticmethod
    def _ink2sum(uG,ni,nj):
        v2nodes = list(uG.neighbors(nj))
        for neigh in uG.neighbors(nj):
            v2nodes+list(uG.neighbors(neigh))
        num = 0
        for node in set(v2nodes):
            num+= uG.degree(node)
        return num
    @staticmethod
    def _k1sum(uG,ni,nj):
        num = 0
        for neigh in uG.neighbors(nj):
            num+= uG.degree(neigh)
        return num
    
    @staticmethod
    def indegree(G,ni,nj,rand_node):
        return G.in_degree(nj),G.in_degree(rand_node)
    @staticmethod
    def outdegree(G,ni,nj,rand_node):
        return G.out_degree(nj),G.out_degree(rand_node)
    @staticmethod
    def degree(G,ni,nj,rand_node):
        return G.degree(nj),G.degree(rand_node)
    @staticmethod
    def continuous_degree(G,ni,nj,rand_node):
        deg = [len(nbrs) for n,nbrs in G.adj.items()]
        ave = np.average(deg)
        return G.degree(nj)/ave,G.degree(rand_node)/ave
    
    @staticmethod
    def ink1sum(g,ni,nj,rand_node):
        #return index_global._ink1sum(g,ni,nj),index_global._ink1sum(g,ni,rand_node)
        return index_global._k1sum(g,ni,nj), index_global._k1sum(g,ni,nj)

    @staticmethod
    def ink2sum(g,ni,nj,rand_node):
        return index_global._ink2sum(g,ni,nj),index_global._ink2sum(g,ni,rand_node)
    
    @staticmethod #eigenvector
    def katz(uG,ni,nj,rand_node):
        katz = nx.katz_centrality_numpy(uG)
        ave = np.average(list(katz.values()))
        return katz[nj]/ave,katz[rand_node]/ave
    @staticmethod #eigenvector
    def eigenvector(uG,ni,nj,rand_node):
        eigvec = nx.eigenvector_centrality_numpy(uG)
        ave = np.average(list(eigvec.values()))
        return eigvec[nj]/ave,eigvec[rand_node]/ave
    @staticmethod #closeness
    def closeness(uG,ni,nj,rand_node):
        a = nx.closeness_centrality(uG,u=nj)
        b = nx.closeness_centrality(uG,u=rand_node)
        return a,b
    @staticmethod
    def information(uG,ni,nj,rand_node): #information centrality
        d = nx.centrality.information_centrality(uG)
        # ave is ok
        return d[nj],d[rand_node]
    @staticmethod
    def betweeness_centrality(uG,ni,nj,rand_node):
        d = nx.betweenness_centrality(uG)
        return d[nj],d[rand_node]
    @staticmethod
    def current_flow_betweenness(uG,ni,nj,rand_node):
        d = nx.current_flow_betweenness_centrality(uG)
        return d[nj],d[rand_node]
    @staticmethod
    def communicability_betweenness(uG,ni,nj,rand_node):
        d = nx.communicability_betweenness_centrality(uG)
        return d[nj],d[rand_node]
    @staticmethod
    def communicability_exp(uG,ni,nj,rand_node):
        c = nx.communicability_exp(uG)
        return c[ni][nj],c[ni][rand_node]
        
        
    @staticmethod
    def harmonic(uG,ni,nj,rand_node):
        d = nx.harmonic_centrality(uG,[nj,rand_node]) 
        return d[nj],d[rand_node]
    
    @staticmethod #link analysis
    def pagerank(g,ni,nj,rand_node):
        s = g.vs.select(name = ni)[0].index
        t = g.vs.select(name = nj)[0].index
        d = g.pagerank(vertices=[s,t])
        return d[0],d[1]
    @staticmethod
    def h_index(uG,ni,nj,rand_node):
        nj_list = []
        rand_node_list = []
        for neigh in uG.neighbors(nj):
            nj_list.append(uG.degree(neigh))
        for neigh in uG.neighbors(rand_node):
            rand_node_list.append(uG.degree(neigh))
        a = index_global._h_index(nj_list)
        b = index_global._h_index(rand_node_list)
        return a,b
    @staticmethod
    def kshell(uG,ni,nj,rand_node):
        # core number
        kshell = nx.core_number(uG)
        # ave = np.average(list(kshell.values()))
        return kshell[nj],kshell[rand_node]
    @staticmethod
    def hits(uG,ni,nj,rand_node):
        d = nx.hits(uG,tol=0.1)[0]
        return d[nj],d[rand_node]


    
class EdgeList(object):
    def __init__(self):
        self.edge_list = []
        self.time_list = []
        self.name = "none"
        self.records = {}
        self.records_random = {}
        # self.functions = {"degree":indegree,"distance":index_local.cutoff_distance}
        self.functions = None 
        # a functions dictionary example is
        # {"degree":indegree,"distance":index_local.cutoff_distance,"distance":index_local.shortest_path_length}
        # ,"katz":index_global.katz, "k1sum":index_global.ink1sum}
        self.un_functions = None
        self.ig_functions = None
        
    def load_index_for_recorded(self,indices_str):
        # index_list = ["degree","distance",distance_cut,]
        loca = index_local()
        all_index = {"indegree":index_global.indegree,
                     "degree":index_global.degree,
                     "continuous_degree":index_global.continuous_degree,
                     "k1sum_uG":index_global.ink1sum,
                     "k2sum":index_global.ink2sum,
                     "katz_uG":index_global.katz,
                     "eigenvector_uG":index_global.eigenvector,
                     "closeness_uG":index_global.closeness,
                     "information_uG":index_global.information,
                     "betweeness_centrality_uG":index_global.betweeness_centrality,
                     "current_flow_betweenness_uG":index_global.current_flow_betweenness,
                     "communicability_betweenness_uG":index_global.communicability_betweenness,
                     "communicability_exp_uG":index_global.communicability_exp,
                     "harmonic_uG":index_global.harmonic,
                     "pagerank_ig":index_global.pagerank,
                     "h_index_uG":index_global.h_index,
                     "kshell_uG":index_global.kshell,
                     "hits_uG": index_global.hits,
                     
                     
                     "shortest_path_length_uG":index_local.shortest_path_length,
                     "distance_cutoff_uG":index_local.cutoff_distance,
                     "resource_allocation_uG":index_local.resource_allocation,
                     "jaccard_coefficient_uG":index_local.jaccard_coefficient,
                     "adamic_adar_index_uG":index_local.adamic_adar_index,
                     "preferential_attachment_uG":index_local.preferential_attachment,
                     "cn_soundarajan_hopcroft_uG":index_local.cn_soundarajan_hopcroft,
                     "ra_index_soundarajan_hopcroft_uG":index_local.ra_index_soundarajan_hopcroft,
                     "within_inter_cluster_uG":index_local.within_inter_cluster,
                     "common_neighbor_centrality_uG":index_local.common_neighbor_centrality,
                     "salton_uG":index_local.salton,
                     "sorensen_uG":index_local.sorensen,
                     "HPI_uG":index_local.HPI,
                     "HDI_uG":index_local.HDI,
                     "LHN1_uG":index_local.LHN1,
                     "average_commute_time_ig":loca.average_commute_time,
                     "cos_ig":loca.cos,
                     "random_walk_with_restart_ig":index_local.random_walk_with_restart,
                     "matrix_forest_index_ig":index_local.matrix_forest_index,
                     "efficiency_uG":index_local.efficiency
                     }
        self.functions = {}
        for s in indices_str:
            self.functions[s] = all_index[s]
            
    def read_edgelist_file(self,file_name):
        # return self.__read_soc_edgelist_file(file_name)
        return self.__read_mathoverflow_file(file_name)
    
    def __read_mathoverflow_file(self,file_name):
        self.name = file_name.split('/')[-1].split(".")[-2]
        fp = open(file_name, "r")
        original_edge_list = []
        erase_num = 0
        while 1:
            line = fp.readline()
            if not line:
                break
            temp = line.split()
            # new_time = int(temp[2])
            
            n0 = int(temp[0])
            n1 = int(temp[1])
            t = int(temp[2])
            
            if n0!= n1:
                original_edge_list.append([n0,n1,t])
            else:
                erase_num +=1
        fp.close()
        legal_edge_num = len(original_edge_list)
        erase_p = float(erase_num)/(erase_num + legal_edge_num)
        self.edge_list = sorted(original_edge_list, key=lambda x:x[2])
        print('in file load, we erase %f' %(erase_p) )
        return self.edge_list

        
    def __read_soc_edgelist_file(self,file_name):
        self.name = 'soc_alpha'
        fp = open(file_name, "r")
        original_edge_list = []
        erase_num = 0
        while 1:
            line = fp.readline()
            if not line:
                break
            temp = line.split(',')
            new_time = int(temp[2])
            
            n0 = int(temp[0])
            n1 = int(temp[1])
            t = int(temp[3])
            
            if n0!= n1:
                original_edge_list.append([n0,n1,t])
            else:
                erase_num +=1
        fp.close()
        legal_edge_num = len(original_edge_list)
        erase_p = float(erase_num)/(erase_num + legal_edge_num)
        self.edge_list = sorted(original_edge_list, key=lambda x:x[2])
        print('in file load, we erase %f' %(erase_p) )
        return self.edge_list
    
    def read_HK(self,alpha,p,N):
        m = 5
        g = nx.complete_graph(m)
        g.to_directed()
        self.name = "HK_network"
        t =1
        
        for node in range(m-1,N):
            weights = np.power([d for n,d in g.degree()],alpha)
            sum_weights = np.sum(weights)
            prob = [float(weight)/sum_weights for weight in weights] 
            c_node = np.random.choice(g.nodes,p=prob)
            rand_node = np.random.choice(g.nodes)
            obvious_BA_node = c_node
            
            g.add_edge(node,c_node) # first node add by preferential degree attachment
            self.edge_list.append((node,c_node,t))
            
            for i in range(1,m,1):
                if_TF = np.random.choice([0,1], p =[1-p,p])
                if if_TF == 1:
                    choose_nodes =list(g.neighbors(obvious_BA_node))
                    for i in g.neighbors(node):
                        if i in choose_nodes:
                            choose_nodes.remove(i)
                    if len(choose_nodes) ==0:
                        degrees = np.array([d for n,d in g.degree()])
                        weights = np.power(degrees ,alpha)
                        sum_weights = np.sum(weights)
                        prob = [float(i) / sum_degree for i in weights]
                        # prob /= prob.sum()
                        c_node = np.random.choice(g.nodes,1,p=prob)
                        obvious_BA_node = c_node
                        g.add_edge(node, c_node)
                        self.edge_list.append((node, c_node,t))
                        t +=1
                    else:  # reall TF step
                        TF_node = np.random.choice(choose_nodes)
                        g.add_edge(node, TF_node)
                        self.edge_list.append((node, TF_node,t)) 
                        t+=1              
                else:
                    degrees = np.array([d for n,d in g.degree()])
                    weights = np.power(degrees, alpha)
                    sum_degree = np.sum(weights)
                    prob = [float(i) / sum_degree for i in weights]
                    c_node = np.random.choice(g.nodes, p=prob)
                    obvious_BA_node = c_node
                    g.add_edge(node, c_node)
                    self.edge_list.append((node, c_node,t))        
                    t+=1        
        pass
        
    def read_BA(self, alpha, N):
        m = 5
        # N  = 5000
        g = nx.complete_graph(m)
        g = g.to_directed()
        #print(g.nodes)
        self.name = "BA_network"
        t = 1
        for node in range(4,N):
            weights = np.power([d for n,d in g.degree()],alpha)
            #print(weights)
            sum_weights = np.sum(weights)
            prob = [float(weight)/sum_weights for weight in weights]
            nodes_chosen = np.random.choice(g.nodes,m,p=prob)
            for node_chosen in nodes_chosen:
                g.add_edge(node, node_chosen)
                self.edge_list.append((node, node_chosen, t))
                t += 1  
        #print(random.choice(list(g.nodes)))
        return g
    
    def evolve_pool(self, max_edge_number = float("inf")):
        G = nx.DiGraph()
        uG = nx.Graph()
        g = igraph.Graph() # three kinds of graph
        print(type(G))
        records = {}
        #records["time"] = []
        self.time_list = []
        records_random ={}
        records_time = {}
        N = min(len(self.edge_list),max_edge_number)

        for name_func in self.functions:
            records[name_func] = []
            records_random[name_func] = []
            records_time[name_func] =0
        records_time['community'] = 0
        print("EdgeList of "+ self.name + " has "+ str(N) +' edges')
        
        # do not need origin network
        # for i in range(200):
        #     edge = self.edge_list[i]
        #     G.add_edge(edge[0],edge[1])
        #network evolve for a period
        num_edges_i = 0
        num_edges_ii = 0
        num_edges_iii = 0
        num_edges_iv = 0
            
        
        pool = mp.Pool()
        
        for i in tqdm(range(N)):
        #for i in range(2*N):
            
            # print("we are recording the edge number ", i)
            edge = self.edge_list[i]
            ni=edge[0]
            nj=edge[1]
            t = edge[2]
            
            if ni in G.nodes() and nj in G.nodes():
                rand_node = np.random.choice([x for x in G.nodes if x != ni]) #random node but not ni
                t0= time.time()
                p_louvain = g.community_multilevel()
                cluster_index = 0
                for cluster in p_louvain:
                    cluster_index += 1
                    for node in cluster:
                        uG.nodes[g.vs[node]["name"]]["community"] = cluster_index
                records_time['community'] += time.time()-t0
                
                num_edges  =  num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                
                for name_func in self.functions: # record the directed graph indices
                    function = self.functions[name_func]
                    t1 = time.time()
                    if name_func[-2:] == "uG":
                        #target_value,rand_value = function(uG,edge[0],edge[1],rand_node)
                        pool.apply_async(function,(uG,edge[0],edge[1],rand_node))
                    elif name_func[-2:] == "ig":
                        pool.apply_async(function,(g,edge[0],edge[1],rand_node)) 
                    else:
                        pool.apply_async(function,(G,edge[0],edge[1],rand_node)) 
                    
                 
                    target_value, rand_value = 0,1
                    records[name_func].append(target_value)
                    records_random[name_func].append(rand_value)
                    t2 = time.time()
                    records_time[name_func] += t2-t1 
                self.records["time"] = num_edges - uG.nodes[nj]["time"]
                self.records_random["time"] = num_edges - uG.nodes[rand_node]["time"]
                self.time_list.append(t) #record the indices and edge
                
                G.add_edge(ni,nj) # build the new edge
                uG.add_edge(ni,nj)
                s = g.vs.select(name = ni)[0].index
                d = g.vs.select(name = nj)[0].index
                g.add_edges([(s,d)])
                num_edges_i +=1 #this is a type i edge, which is the edges that we consider most important 
                
            elif ni not in G.nodes and nj in G.nodes:
                
                G.add_edge(ni,nj) # ni is a new node, but nj is not a new node
                uG.add_edge(ni,nj)
                g.add_vertices([ni])
                s = g.vs.select(name = ni)[0].index
                d = g.vs.select(name = nj)[0].index
                g.add_edges([(s,d)])
                uG.nodes[ni]["time"] = num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                num_edges_ii +=1
            elif ni in G.nodes and nj not in G.nodes:
            
                G.add_edge(ni,nj)
                uG.add_edge(ni,nj)
                g.add_vertices([nj])
                s = g.vs.select(name = ni)[0].index
                d = g.vs.select(name = nj)[0].index
                g.add_edges([(s,d)])
                uG.nodes[nj]["time"] = num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                num_edges_iii +=1 # ni is an old node, but nj is a new node
            else:
                G.add_edge(ni,nj)
                uG.add_edge(ni,nj)
                g.add_vertices([ni,nj])
                s = g.vs.select(name = ni)[0].index
                d = g.vs.select(name = nj)[0].index
                g.add_edges([(s,d)])
                uG.nodes[ni]["time"] = num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                uG.nodes[nj]["time"] = num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                num_edges_iv +=1 # both ni, nj are new node
            if num_edges_i >N :
                break
        pool.close()
        pool.join()
        self.records = records
        self.records_random = records_random
        #print(num_edges_i,num_edges_ii,num_edges_iii,num_edges_iv)
        print("type i edges is " + str(num_edges_i/N) +" in record edges!")
        print("type ii edges is " + str(num_edges_ii/N) +" in record edges!") 
        print("type iii edges is " + str(num_edges_iii/N) +" in record edges!")
        print("type iv edges is " + str(num_edges_iv/N) +" in record edges!")
        
        print(sorted(records_time.items(), key = lambda kv:(kv[1], kv[0]))) #打印所有方法使用时间
        return G
      
    def evolve(self, max_edge_number = float("inf")):
        G = nx.DiGraph()
        uG = nx.Graph()
        g = igraph.Graph() # three kinds of graph
        print(type(G))
        records = {}
        #records["time"] = []
        self.time_list = []
        records_random ={}
        records_time = {}
        N = min(len(self.edge_list),max_edge_number)

        for name_func in self.functions:
            records[name_func] = []
            records_random[name_func] = []
            records_time[name_func] =0
        records_time['community'] = 0
        print("EdgeList of "+ self.name + " has "+ str(N) +' edges')
        
        # do not need origin network
        # for i in range(200):
        #     edge = self.edge_list[i]
        #     G.add_edge(edge[0],edge[1])
        #network evolve for a period
        num_edges_i = 0
        num_edges_ii = 0
        num_edges_iii = 0
        num_edges_iv = 0
            
        
        
        for i in tqdm(range(N)):
        #for i in range(2*N):
            
            # print("we are recording the edge number ", i)
            edge = self.edge_list[i]
            ni=edge[0]
            nj=edge[1]
            t = edge[2]
            
            if ni in G.nodes() and nj in G.nodes():
                rand_node = np.random.choice([x for x in G.nodes if x != ni]) #random node but not ni
                t0= time.time()
                p_louvain = g.community_multilevel()
                cluster_index = 0
                for cluster in p_louvain:
                    cluster_index += 1
                    for node in cluster:
                        uG.nodes[g.vs[node]["name"]]["community"] = cluster_index
                records_time['community'] += time.time()-t0
                
                num_edges  =  num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                
                for name_func in self.functions: # record the directed graph indices
                    function = self.functions[name_func]
                    t1 = time.time()
                    if name_func[-2:] == "uG":
                        target_value,rand_value = function(uG,edge[0],edge[1],rand_node)
                        records[name_func].append(target_value)
                        records_random[name_func].append(rand_value)
                    elif name_func[-2:] == "ig":
                        target_value,rand_value = function(g,edge[0],edge[1],rand_node)
                        records[name_func].append(target_value)
                        records_random[name_func].append(rand_value)
                    else:
                        target_value,rand_value = function(G,edge[0],edge[1],rand_node)
                        records[name_func].append(target_value)
                        records_random[name_func].append(rand_value)
                    t2 = time.time()
                    records_time[name_func] += t2-t1 
                self.records["time"] = num_edges - uG.nodes[nj]["time"]
                self.records_random["time"] = num_edges - uG.nodes[rand_node]["time"]
                self.time_list.append(t) #record the indices and edge
                
                G.add_edge(ni,nj) # build the new edge
                uG.add_edge(ni,nj)
                s = g.vs.select(name = ni)[0].index
                d = g.vs.select(name = nj)[0].index
                g.add_edges([(s,d)])
                num_edges_i +=1 #this is a type i edge, which is the edges that we consider most important 
                
            elif ni not in G.nodes and nj in G.nodes:
                
                G.add_edge(ni,nj) # ni is a new node, but nj is not a new node
                uG.add_edge(ni,nj)
                g.add_vertices([ni])
                s = g.vs.select(name = ni)[0].index
                d = g.vs.select(name = nj)[0].index
                g.add_edges([(s,d)])
                uG.nodes[ni]["time"] = num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                num_edges_ii +=1
            elif ni in G.nodes and nj not in G.nodes:
            
                G.add_edge(ni,nj)
                uG.add_edge(ni,nj)
                g.add_vertices([nj])
                s = g.vs.select(name = ni)[0].index
                d = g.vs.select(name = nj)[0].index
                g.add_edges([(s,d)])
                uG.nodes[nj]["time"] = num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                num_edges_iii +=1 # ni is an old node, but nj is a new node
            else:
                G.add_edge(ni,nj)
                uG.add_edge(ni,nj)
                g.add_vertices([ni,nj])
                s = g.vs.select(name = ni)[0].index
                d = g.vs.select(name = nj)[0].index
                g.add_edges([(s,d)])
                uG.nodes[ni]["time"] = num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                uG.nodes[nj]["time"] = num_edges_i + num_edges_ii + num_edges_iii + num_edges_iv + 1
                num_edges_iv +=1 # both ni, nj are new node
            if num_edges_i >N :
                break
        self.records = records
        self.records_random = records_random
        #print(num_edges_i,num_edges_ii,num_edges_iii,num_edges_iv)
        print("type i edges is " + str(num_edges_i/N) +" in record edges!")
        print("type ii edges is " + str(num_edges_ii/N) +" in record edges!") 
        print("type iii edges is " + str(num_edges_iii/N) +" in record edges!")
        print("type iv edges is " + str(num_edges_iv/N) +" in record edges!")
        
        print(sorted(records_time.items(), key = lambda kv:(kv[1], kv[0]))) #打印所有方法使用时间
        return G
    
    def output_records(self):
        data = pd.DataFrame.from_dict(self.records)
        data_random = pd.DataFrame.from_dict(self.records_random)
        data_time = pd.DataFrame.from_dict({"time":self.time_list})
        # data.to_csv(self.name+"indeces")
        # data_random.to_csv(self.name+"random_indeces")
        # data_time.to_csv(self.name+"time")
        writer = pd.ExcelWriter("./experimental_results/"+self.name+'_all.xlsx')
        data.to_excel(writer,sheet_name = "real",index=False)
        data_random.to_excel(writer,sheet_name = "random",index=False)
        data_time.to_excel(writer,sheet_name ="time",index=False)
        writer.save()
        writer.close()
        # result = pd.concat([data,data_random,data_time], keys=['real', 'random', 'time'])
        # result.to_csv(self.name +"all.csv")
        # self.records
        pass
    
    def load_records(self,file_name):
        data_real = pd.read_excel(file_name,sheet_name = "real")
        self.records = {k:list(v.values()) for k,v in data_real.to_dict().items()}
        data_real = pd.read_excel(file_name,sheet_name = "random")
        self.records_random = {k:list(v.values()) for k,v in data_real.to_dict().items()}
        data_real = pd.read_excel(file_name,sheet_name = "time")
        self.time_list = {k:list(v.values()) for k,v in data_real.to_dict().items()}  
        
    def cut_records(self,normalize_indices,smooth_length):
        # just cut off 
        for index_name in normalize_indices:
            self.records[index_name] = (np.array(self.records[index_name]))[smooth_length-1:-smooth_length+1] #normalize and cut_off
            self.records_random[index_name] = (np.array(self.records_random[index_name]))[smooth_length-1:-smooth_length+1]  #add to avoid zero
            #self.time_list = np.array(self.time_list)[smooth_length-1:-smooth_length+1]
    
    def add_small_number_avoid_divide_by_zero(self,avoid_zero_indices,epsilon):
        for index_name in avoid_zero_indices:
            self.records[index_name] = np.array(self.records[index_name])+epsilon
            self.records_random[index_name] = np.array(self.records_random[index_name]) +epsilon
            
        
    def cut_and_smooth_normalize_records(self,normalize_indices,smooth_length):
        # to check
        # do three things: smooth, normalize and cutoff
        weights = np.ones(smooth_length)/smooth_length
        for index_name in normalize_indices:
            x = np.array(self.records_random[index_name]) 
            x_ave = np.convolve(x,weights,'same') #smooth
            
            self.records[index_name] = ((np.array(self.records[index_name]))/x_ave)[smooth_length-1:-smooth_length+1] #normalize and cut_off and add 1 to avoid zero
            self.records_random[index_name] = (np.array(self.records_random[index_name])/x_ave)[smooth_length-1:-smooth_length+1] 
            #self.time_list = np.array(self.time_list)[smooth_length-1:-smooth_length+1]
    

def merge(indices,indices_random):
    index_merge = []
    if_random = []
    for i in range(len(indices)):
        index_merge.append(indices[i])
        index_merge.append(indices_random[i])
        if_random.append(0)
        if_random.append(1)
        
    return np.array(index_merge), np.array(if_random)
                                           
def first_plot():
    el = EdgeList()
    file_name = './BA_network_all.xlsx'
    el.load_records(file_name)
    #el.smooth_and_normalize_records(sl_normalize_indices=['degree'],smooth_length=100)
    degree,if_rand = merge(el.records["degree"],el.records_random["degree"])
    distances, if_rand = merge(el.records["distance"],el.records_random["distance"])
    print(if_rand.dtype=='int64')
    print(degree.dtype=='float64')
    #print(len(distances),len(if_rand),len(degree))
    
    step_num = len(if_rand)
    step_size = step_num // 10
    edge_nums = np.arange(step_size, len(if_rand),step_size)
    ent1 = []
    ent2 = []
    ent3 = []
    ent4 = []
    ent5 = []
    ent6 = []
    ent0 = []
    ent7 = []
    ents = {}
    
    for i in edge_nums:
        if_rand_cut = if_rand[i-step_size:i]
        degree_cut = degree[i-step_size:i]
        distances_cut = distances[i-step_size:i]
        #print(len(if_rand_cut),len(degree_cut))
        
        k = int(np.sqrt(len(degree_cut)))
        slice_num = 20
        
        slice_num = int(np.power(len(degree_cut),1.0/3))
        # print("length of datas, number of spaces",len(degree_cut),slice_num)
        #print(k)
        ent1.append(ep.midcut(if_rand_cut, degree_cut,slice_num=slice_num))
        ent2.append(ep.midc(if_rand_cut, degree_cut,k=20))
        
        # ent0.append(ep.cmidcutd(if_rand_cut, degree_cut,distances_cut))
        k = 20
        a,b = ep.cmiddc(if_rand_cut,distances_cut,degree_cut,k=k)
        ent6.append(a)
        ent7.append(b)

        ent3.append(ep.cmiddcut(if_rand_cut, distances_cut, degree_cut))
        if_rand_cut_copy = if_rand_cut.copy()
        random.shuffle(if_rand_cut_copy)
        
        ent4.append(ep.cmiddcut(if_rand_cut_copy,distances_cut,degree_cut))
        random.shuffle(distances_cut)
        ent5.append(ep.cmiddcut(if_rand_cut,distances_cut,degree_cut))
        
        #ent5.append(ep.cmidcutd(if_rand_cut,degree_cut,distances_cut))
        #print(a.all(),b.all(),c.all())
        #ent5.append(ep.cmicut(if_rand_cut, degree_cut,distances_cut))
        #ent6.append(ep.cmi(if_rand_cut, degree_cut,distances_cut))
    
    print('len')
    print(len(ent6))
    print(len(ent0))
    # print(ent6)

    

    # x_vec = np.array([[s] for s in if_rand ])
    # y_vec = np.array([[s] for s in degree ])

    
    # print("midc",ee.midc(x_vec,y_vec,base=2,k=20))
    # print("midc",ep.midc(if_rand,degree,k=40,base=2))
    # print(ee.shuffle_test(ee.midc,x_vec,y_vec,base=2,k=10))
    
    # print(ep.midc(if_rand,degree,k=20,base=2))
    
    slice_num = int(np.power(len(degree),1.0/2))    
    
    print(ep.midd(if_rand,distances))
    x_vec = np.transpose(np.array([if_rand]))
    y_vec = np.transpose(np.array([distances]))
    z_vec = np.transpose(np.array([ep.slice(degree, slice_num)]))
    print(ee.shuffle_test(ee.midd,x_vec,y_vec,base=2))
    print(ee.shuffle_test(ee.midd,x_vec,z_vec,base=2))
    print(ee.cmidd(x_vec,z_vec,y_vec,base=2))
    # print(ee.shuffle_test(ee.cmidd,x_vec,y_vec,z_vec,base=2))

    # plt.plot(ent0,label = r'0')
    # plt.plot(ent6,label = r"6")
    plt.plot(ent7,label = r"7")
    plt.plot(ent1, 'vb-',label=r'1')
    plt.plot(ent2, 'or--',label=r'2')
    plt.plot(ent3,label=r'3')
    plt.plot(ent4,label=r'4')
    plt.plot(ent5,label=r'5')
    # plt.plot(edge_nums, ent5, 'b',label='5')
    # plt.plot(edge_nums, ent6, 'b--',label='6')
    plt.legend()
    plt.show()

def first_plot2():
    el = EdgeList()
    el.load_records("file_name")
    #el.smooth_and_normalize_records(sl_normalize_indices=['degree'],smooth_length=100)
    degree,if_rand = merge(el.records["degree"],el.records_random["degree"])
    distances, if_rand = merge(el.records["distance"],el.records_random["distance"])
    #print(len(distances),len(if_rand),len(degree))
    
    step_num = len(if_rand)
    step_size = step_num // 10
    edge_nums = np.arange(step_size, len(if_rand),step_size)

    ents = {}
    k_max = int(np.sqrt(step_size))
    ks = range(20,k_max)
    for k in ks:
            ents[k] = []
    ent0 = []
    for i in edge_nums:
        if_rand_cut = if_rand[i-step_size:i]
        degree_cut = degree[i-step_size:i]
        distances_cut = distances[i-step_size:i]
        a,b,c,d = ep.entropy_package(if_rand_cut,degree_cut,distances_cut)
        #print(len(if_rand_cut),len(degree_cut))
        #print(k)
        
        for k in ks:
            ents[k].append(ep.midcut(if_rand_cut, degree_cut,slice_num=k))
            
        ent0.append(ep.midc(if_rand_cut, degree_cut,k=100))
        #ent5.append(ep.cmidcutd(if_rand_cut,degree_cut,distances_cut))
        #print(a.all(),b.all(),c.all())
        #ent5.append(ep.cmicut(if_rand_cut, degree_cut,distances_cut))
        #ent6.append(ep.cmi(if_rand_cut, degree_cut,distances_cut))
    print('len',step_size)
    print('k_max', k_max)
    for k in ks:
        plt.plot(ents[k], label = str(k))

    plt.plot(ent0, 'vb-',label=r'continuous')
    # plt.plot(edge_nums, ent5, 'b',label='5')
    # plt.plot(edge_nums, ent6, 'b--',label='6')
    plt.legend()
    plt.show()

def test():
    el = EdgeList()
    g = el.read_BA(1.0,300)
    #print(el.edge_list)
    # li = el._evolve_test()
    print(g.degree)
    # print(li)
     

def analysis_index(file_name):
    el = EdgeList()
    
    el.load_records(file_name)
    index4normalize = ["indegree",
                       "degree",
                       "continuous_degree",
                       "k1sum",
                       "k2sum",
                       "katz_uG",
                       #running time"eigenvector_uG",
                       "closeness_uG",
                       ## connected graph "information_uG",
                       ## running time "betweeness_centrality_uG",
                       ## connected graph "current_flow_betweenness_uG",
                       ## "communicability_betweenness_uG",
                       "harmonic_uG",
                       #float division by zero "rich_club_coefficient_uG",
                       "pagerank_ig",
                       "h_index_uG",
                       "kshell_uG",
                       
                       #"resource_allocation_uG",
                       #"jaccard_coefficient_uG",
                       #"adamic_adar_index_uG",
                       "preferential_attachment_uG",
                       #"cn_soundarajan_hopcroft_uG",
                       #"ra_index_soundarajan_hopcroft_uG",
                       #"within_inter_cluster_uG",
                       "common_neighbor_centrality_uG",
                       "efficiency_uG"
                       ]
    
    index4cutoff = ['distance_uG','distance_cutoff_uG']
    index4addepsilon = ["resource_allocation_uG",
                    "jaccard_coefficient_uG","adamic_adar_index_uG","cn_soundarajan_hopcroft_uG",
                   "ra_index_soundarajan_hopcroft_uG","within_inter_cluster_uG"]
    
    el.cut_records(normalize_indices=index4cutoff,smooth_length=300)
    el.add_small_number_avoid_divide_by_zero(avoid_zero_indices=index4addepsilon,epsilon=1)    
    el.cut_and_smooth_normalize_records(normalize_indices=index4normalize+index4addepsilon,smooth_length=300)
    

    total_slice_num = 5
    ordered_slice_num = 4 #1,2,3,4
    index_min =0
    index_max =1
    length = len(el.records["degree"])
    
    print(type(el.records["distance_uG"][0]))
    print(type(el.records["distance_cutoff_uG"][1]))
    
    if ordered_slice_num !=4:
        index_min = ordered_slice_num*int(length/total_slice_num)
        index_max = (ordered_slice_num+1)*int(length/total_slice_num)
    else:
        index_min = ordered_slice_num*int(length/total_slice_num)
        index_max = length - 1
        
        
    
    
    index_all_list = index4normalize+index4cutoff
    index_data ={}
    degree,if_rand = merge(el.records["degree"][index_min:index_max],el.records_random["degree"][index_min:index_max])
    
    index_data["if_rand"] = if_rand
    for index in index_all_list:
        index_records, if_rand =merge(el.records[index][index_min:index_max],el.records_random[index][index_min:index_max])
        index_data["if_rand"] = if_rand
        index_data[index] = index_records
        
        
    
    # degree,if_rand = merge(el.records["degree"],el.records_random["degree"])
    # distance, if_rand = merge(el.records["distance"],el.records_random["distance"])
    # indegree, if_rand = merge(el.records["indegree"],el.records_random["indegree"])
    # distance_cutoff, if_rand = merge(el.records["distance_cutoff"],el.records_random["distance_cutoff"])
    # katz, if_rand = merge(el.records["katz"],el.records_random["katz"])
    # k1sum, if_rand = merge(el.records["k1sum"],el.records_random["k1sum"])
    # k2sum, if_rand = merge(el.records["k2sum"],el.records_random["k2sum"])

    
    
    # index_data = {"if_rand":if_rand,"degree":degree,"indegree":indegree,"distance_cutoff":distance_cutoff,
    #               "distance":distance,"katz":katz,"k1sum":k1sum,"k2sum":k2sum}
    # index_nodes = ["degree"]
    # index_all = set(["degree","distance_cutoff","distance","katz","indegree","k1sum","k2sum"])
    index_all = set(index_all_list)
    
    def aggregative_discovery():
        K = set([]) 
        x = 100
        p = None
        
        while x > 0:
            if p != None:
                K.add(p)
            print("K is", K)
            max_x = 0
            max_j = None
            index_check = index_all.difference(K)
            print("index for check is ",index_check)
            for j in index_check:
                v,_,_,_,_,multi = ep.casual_entropy("if_rand",j,K,index_data)
                if v>max_x:
                    max_x = v
                    max_j = j
                else:
                    # print(j+" useful value",v)
                    pass
            x = max_x
            p = max_j
        return K
    
    def aggregative_discovery_by_multi():
        K = set([]) 
        x = 100
        p = None
        
        while x > 0:
            if p != None:
                K.add(p)
            print("K is", K)
            max_multi = 0
            max_j = None
            index_check = index_all.difference(K)
            print("index for check is ",index_check)
            for j in index_check:
                v,_,_,_,_,multi = ep.casual_entropy("if_rand",j,K,index_data)
                if multi>max_multi:
                    max_multi = multi
                    max_j = j
                else:
                    #print(j+" useful value",v)
                    pass
            x = max_multi
            p = max_j
        return K 
    
    def progressive_removal(K):
        for j in K:
            L = K.copy()
            L.difference(set([j]))
            v,_,_,_,_,_ = ep.casual_entropy("if_rand",j,L,index_data)
            if v <=0:
                K.difference(j)
        return K
    
    K  = aggregative_discovery()
    print(file_name+" 's final K is",progressive_removal(K))

    
    sent_notifer()
    # useful_value,v,ave,(ci0,ci1),if_large_zero, multi =ep.casual_entropy('if_rand','distance',set(['indegree','degree']),index_data)
    # print(useful_value,v,ave,(ci0,ci1),if_large_zero,multi)
    pass


def test_casual_entropy():
    el = EdgeList()
    file_name = './BA_network_all.xlsx'
    el.load_records(file_name)
    
def test_numpy_c():
    x = np.array([[1],[2]])
    y = np.array([[3],[4]])
    z = np.array([[5],[6]])
    a = [x,y]
    
    print(np.c_[(x,y,z)])
    print(np.concatenate((x,y,z),axis=1 ))
    pass

def test_all_network():
    pass

if __name__ == "__main__":
    
    file_name = './real_evolving_networks/hepph.txt'
    # t0 =time.time()
    network_evolve(file_name,max_edge_number=300)
    # print("test time ", time.time()-t0)
    # el =EdgeList()
    # el.load_records()
    # file_name = "./experimental_result"
    
    # file_name = "./experimental_results/hepph_all.xlsx"
    # analysis_index(file_name)
    
   
    # first_plot()
    
    # first_plot2()
    # el = EdgeList()
    # g = el.read_BA(1.0,1000)
    
    # print(list(g.predecessors(2)))
    


    
    # print(index_global.katz(g,1,2))
    # for a,b,c in nx.resource_allocation_index(g, [(1,2)]):
    #     print(a,b,c)
    # print(list(nx.resource_allocation_index(g, [(1,2)])))
    # #print(el.edge_list)


    # el.smooth_and_normalize_records(sl_normalize_indices=['degree'],smooth_length=100)
    # print("records",el.records)
    # print("records_random",el.records_random)
    
    #print(nx.degree(dg,50),nx.in_degree(dg,50),nx.out_degree(dg,50))
   # print(el.records)