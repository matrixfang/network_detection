import networkx as nx
import numpy as np
import random 
import matplotlib.pyplot as plt
import entropy as ep
import pandas as pd
from npeet import entropy_estimators as ee
import os
from tqdm import tqdm

def sent_notifer():
    command = "terminal-notifier -activate 'com.microsoft.VSCode' -sender 'com.microsoft.VSCode' -title 'python notifer' -message 'Your python program has finished!' "
    os.system(command)
    os.system('say "your program has finished"')
    # com.apple.Safari
    # 
    pass

    


class index_local(object):
    def __init__(self):
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
    def cutoff_distance(g,ni,nj,rand_node):
        return index_local._cutoff_distance(g,ni,nj),index_local._cutoff_distance(g,ni,rand_node)
        
    @staticmethod
    def _shortest_path_length(g,ni,nj):
        try:
            length = nx.shortest_path_length(g, source=ni, target=nj)
        except nx.exception.NetworkXNoPath:
            return 100
        else:
            return length
    @staticmethod
    def shortest_path_length(g,ni,nj,rand_node):
        return index_local._shortest_path_length(g,ni,nj),index_local._shortest_path_length(g,ni,rand_node)
        
    # bascially all the link prediction index is locally index, can be viewed as some kind of similarity
    

class index_global(object):
    def __init__(self):
        # ni does not affect the value of index_global
        pass
    @staticmethod
    def indegree(g,ni,nj,rand_node):
        return g.in_degree(nj),g.in_degree(rand_node)
    @staticmethod
    def outdegree(g,ni,nj,rand_node):
        return g.out_degree(nj),g.out_degree(rand_node)
    @staticmethod
    def degree(g,ni,nj,rand_node):
        return g.degree(nj),g.degree(rand_node)
    @staticmethod
    def continuous_degree(g,ni,nj,rand_node):
        ave = np.average(list(g.degree().values()))
        return g.degree(nj)/ave,g.degree(rand_node)/ave
    @staticmethod
    def katz(g,ni,nj,rand_node):
        katz = nx.katz_centrality_numpy(g)
        ave = np.average(list(katz.values()))
        return katz[nj]/ave,katz[rand_node]/ave
    @staticmethod
    def eigenvector(g,ni,nj,rand_node):
        eigvec = nx.eigenvector_centrality_numpy(g)
        ave = np.average(list(eigvec.values()))
        return eigvec[nj]/ave,eigvec[rand_node]/ave
    
    @staticmethod
    def kshell(g,ni,nj,rand_node):
        # core number
        kshell = nx.core_number(g)
        # ave = np.average(list(kshell.values()))
        return kshell[nj],kshell[rand_node]
    @staticmethod
    def _ink1sum(g,ni,nj):
        num = 0
        for neigh in g.predecessors(nj):
            num += g.in_degree(neigh)
        return num
    @staticmethod
    def ink1sum(g,ni,nj,rand_node):
        return index_global._ink1sum(g,ni,nj),index_global._ink1sum(g,ni,rand_node)
    @staticmethod
    def _ink2sum(g,ni,nj):
        v2nodes = list(g.predecessors(nj))
        for neigh in g.predecessors(nj):
            v2nodes+list(g.predecessors(neigh))
        num = 0
        for node in set(v2nodes):
            num+= g.in_degree(node)
        return num
    @staticmethod
    def ink2sum(g,ni,nj,rand_node):
        return index_global._ink2sum(g,ni,nj),index_global._ink2sum(g,ni,rand_node)

    
class EdgeList(object):
    def __init__(self):
        self.edge_list = []
        self.time_list = []
        self.name = "none"
        self.records = {}
        self.records_random = {}
        # self.functions = {"degree":indegree,"distance":index_local.cutoff_distance}
        self.functions = None 
        # {"degree":indegree,"distance":index_local.cutoff_distance,"distance":index_local.shortest_path_length}
        # ,"katz":index_global.katz, "k1sum":index_global.ink1sum}
        
    def load_index_for_recorded(self,indices_str):
        # index_list = ["degree","distance",distance_cut,]
        all_index = {"indegree":index_global.indegree,"degree":index_global.degree,"distance_cutoff":index_local.cutoff_distance,"distance":index_local.shortest_path_length,
                        "katz":index_global.katz,"k1sum":index_global.ink1sum,"k2sum":index_global.ink2sum}
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
    
    def _evolve_test(self, max_edge_number = float("inf")):
        g = nx.DiGraph()
        print(type(g))
        records = {}
        #records["time"] = []
        self.time_list = []
        records_random ={}
        N = min(len(self.edge_list),max_edge_number)
        print("EdgeList of "+ self.name + " has "+ str(N) +' edges')
        
        li = []
        
        for i in range(200):
            edge = self.edge_list[i]
            g.add_edge(edge[0],edge[1])
        #network evolve for a period
        num_edges_i = 0
        num_edges_ii = 0
        num_edges_iii = 0
        num_edges_iv = 0
        for i in tqdm(range(200, N)):
            # print("we are recording the edge number ", i)
            edge = self.edge_list[i]
            ni=edge[0]
            nj=edge[1]
            t = edge[2]
            if ni in g.nodes() and nj in g.nodes():
                rand_node = np.random.choice(g.nodes)
                a,b = index_global.katz(g,edge[0],edge[1],rand_node)
                li.append(b)
                # for name_func in self.functions:
                #     function = self.functions[name_func]
                #     
                #     records[name_func].append(function(g,edge[0],edge[1]))
                #     records_random[name_func].append(function(g,edge[0],rand_node))
                # self.time_list.append(t)
                g.add_edge(ni,nj) # record and then record the indeces
                num_edges_i +=1 #this is a type i edge, which is the edges that we consider most important 
            elif ni not in g.nodes and nj in g.nodes:
                g.add_edge(ni,nj) # ni is a new node, but nj is not a new node
                num_edges_ii +=1
            elif ni in g.nodes and nj not in g.nodes:
                g.add_edge(ni,nj)
                num_edges_iii +=1 # ni is an old node, but nj is a new node
            else:
                g.add_edge(ni,nj)
                num_edges_iv +=1 # both ni, nj are new node
        

        #print(num_edges_i,num_edges_ii,num_edges_iii,num_edges_iv)
        print("type i edges is " + str(float(num_edges_i)/N) +" in record edges!")
        return li
    
    def evolve(self, max_edge_number = float("inf")):
        g = nx.DiGraph()
        print(type(g))
        records = {}
        #records["time"] = []
        self.time_list = []
        records_random ={}
        N = min(len(self.edge_list),max_edge_number)
        for name_func in self.functions:
            records[name_func] = []
            records_random[name_func] = []
        print("EdgeList of "+ self.name + " has "+ str(N) +' edges')
        
        for i in range(200):
            edge = self.edge_list[i]
            g.add_edge(edge[0],edge[1])
        #network evolve for a period
        num_edges_i = 0
        num_edges_ii = 0
        num_edges_iii = 0
        num_edges_iv = 0
        for i in tqdm(range(200, N)):
            # print("we are recording the edge number ", i)
            edge = self.edge_list[i]
            ni=edge[0]
            nj=edge[1]
            t = edge[2]
            if ni in g.nodes() and nj in g.nodes():
                rand_node = np.random.choice(g.nodes)
                for name_func in self.functions:
                    function = self.functions[name_func]
                    function(g,edge[0],edge[1],rand_node)
                    target_value,rand_value = function(g,edge[0],edge[1],rand_node)
                    records[name_func].append(target_value)
                    records_random[name_func].append(rand_value)
                self.time_list.append(t)
                g.add_edge(ni,nj) # record and then record the indeces
                num_edges_i +=1 #this is a type i edge, which is the edges that we consider most important 
            elif ni not in g.nodes and nj in g.nodes:
                g.add_edge(ni,nj) # ni is a new node, but nj is not a new node
                num_edges_ii +=1
            elif ni in g.nodes and nj not in g.nodes:
                g.add_edge(ni,nj)
                num_edges_iii +=1 # ni is an old node, but nj is a new node
            else:
                g.add_edge(ni,nj)
                num_edges_iv +=1 # both ni, nj are new node
        
        self.records = records
        self.records_random = records_random
        #print(num_edges_i,num_edges_ii,num_edges_iii,num_edges_iv)
        print("type i edges is " + str(float(num_edges_i)/N) +" in record edges!")
        return g
    
    def output_records(self):
        data = pd.DataFrame.from_dict(self.records)
        data_random = pd.DataFrame.from_dict(self.records_random)
        data_time = pd.DataFrame.from_dict({"time":self.time_list})
        # data.to_csv(self.name+"indeces")
        # data_random.to_csv(self.name+"random_indeces")
        # data_time.to_csv(self.name+"time")
        writer = pd.ExcelWriter("./"+self.name+'_all.xlsx')
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
        
    def cutoff_records(self,normalize_indices,smooth_length):
        # just cut off 
        for index_name in normalize_indices:
            self.records[index_name] = (np.array(self.records[index_name]))[smooth_length-1:-smooth_length+1] #normalize and cut_off
            self.records_random[index_name] = (np.array(self.records_random[index_name]))[smooth_length-1:-smooth_length+1]
            #self.time_list = np.array(self.time_list)[smooth_length-1:-smooth_length+1]
        
    def smooth_and_normalize_records(self,normalize_indices,smooth_length): #to check
        # do three things: smooth, normalize and cutoff
        weights = np.ones(smooth_length)/smooth_length
        for index_name in normalize_indices:
            x = self.records_random[index_name]
            x_ave = np.convolve(x,weights,'same') #smooth
            self.records[index_name] = (np.array(self.records[index_name])/x_ave)[smooth_length-1:-smooth_length+1] #normalize and cut_off
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
     
def network_evolve(file_name):
    el = EdgeList() #build edgelist object
    index4normalize = ["indegree","degree","katz","k1sum","k2sum"]
    index4cutoff = ['distance','distance_cutoff']
    el.load_index_for_recorded(index4normalize+index4cutoff)
    # g = el.read_BA(1.0,500)
    # g = el.read_HK(1.0,0.05,1000)
    # g = el.read_HK(1.0,0.05,1000)
    
    el.read_edgelist_file(file_name)
    #print(el.edge_list)
    dg = el.evolve(max_edge_number = 5000)
    el.smooth_and_normalize_records(normalize_indices=index4normalize,smooth_length=100)
    el.cutoff_records(normalize_indices=index4cutoff,smooth_length=100)
    print(el.records)
    print(el.records_random)
    el.output_records()
    sent_notifer()
    # el.load_records("BA_network_all.xlsx")

def analysis_index(file_name):
    el = EdgeList()
    
    el.load_records(file_name)
    degree,if_rand = merge(el.records["degree"],el.records_random["degree"])
    distance, if_rand = merge(el.records["distance"],el.records_random["distance"])
    indegree, if_rand = merge(el.records["indegree"],el.records_random["indegree"])
    distance_cutoff, if_rand = merge(el.records["distance_cutoff"],el.records_random["distance_cutoff"])
    katz, if_rand = merge(el.records["katz"],el.records_random["katz"])
    k1sum, if_rand = merge(el.records["k1sum"],el.records_random["k1sum"])
    k2sum, if_rand = merge(el.records["k2sum"],el.records_random["k2sum"])
    
    
    index_data = {"if_rand":if_rand,"degree":degree,"indegree":indegree,"distance_cutoff":distance_cutoff,
                  "distance":distance,"katz":katz,"k1sum":k1sum,"k2sum":k2sum}
    index_nodes = ["degree"]
    index_all = set(["degree","distance_cutoff","distance","katz","indegree","k1sum","k2sum"])
    
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

if __name__ == "__main__":
    
    file_name = './real_evolving_networks/sx-mathoverflow-a2q.txt'
    network_evolve(file_name)
    # file_name = './hepth_all.xlsx'
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