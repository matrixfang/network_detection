from evolving_mechanism import *

def network_evolve(file_name,max_edge_number=float('inf')):
    el = EdgeList() #build edgelist object
    index4normalize = ["indegree",
                       "degree",
                       "continuous_degree",
                       "k1sum_uG",
                       "k2sum",
                       "katz_uG",
                       "eigenvector_ig",
                       "closeness_uG",
                       ## connected graph "information_uG",
                       ## running time "betweeness_centrality_uG",
                       ## connected graph "current_flow_betweenness_uG",
                       #"communicability_betweenness_uG",
                       #time "communicability_exp_uG",
                       
                       "harmonic_uG",
                       "pagerank_ig",
                       "h_index_uG",
                       "kshell_uG",
                       "hub_score_ig",
                       "add_time_uG",
                       #"hits_uG",
                       
                       
                       #"shortest_path_length_uG",
                       "resource_allocation_uG",
                       "jaccard_coefficient_uG",
                       "adamic_adar_index_uG",
                       
                       "cn_soundarajan_hopcroft_uG",
                       "ra_index_soundarajan_hopcroft_uG",
                       "within_inter_cluster_uG",
                       
                       "preferential_attachment_uG",
                       "common_neighbor_centrality_uG",
                       "salton_uG",
                       "sorensen_uG",
                       "HPI_uG",
                       "HDI_uG",
                       "LHN1_uG",
                       #"average_commute_time_ig",
                       #time"cos_ig",
                       "random_walk_with_restart_ig",
                       #"matrix_forest_index_ig",
                       "efficiency_uG"
                       ]
    
    index4cutoff = ["shortest_path_length_uG",'distance_cutoff_uG']
    el.load_index_for_recorded(index4normalize+index4cutoff)
    # g = el.read_BA(1.0,500)
    # g = el.read_HK(1.0,0.05,1000)
    # g = el.read_HK(1.0,0.05,1000)
    
    el.read_edgelist_file(file_name)
    print("total degelist number ", len(el.edge_list))
    
    #print(el.edge_list)
    dg = el.evolve(max_edge_number = max_edge_number)
    smooth_length = 100 #cut records at begin and end with 100
    el.output_records()
    
    return el
    
    
     
    # el.cut_records(normalize_indices=index4cutoff,smooth_length=50)
    # el.cut_and_smooth_normalize_records(normalize_indices=index4normalize,smooth_length=50)

    # sent_notifer()
    # el.load_records("BA_network_all.xlsx")


def after_analysis_index(file_name):
    el = EdgeList()
    
    el.load_records(file_name)
    index4normalize = ["indegree",
                       "degree", 
                       "k1sum_uG",

                       #float division by zero "rich_club_coefficient_uG",
                       "pagerank_ig",
                       
                       "HDI_uG",
                       "LHN1_uG",
                       "random_walk_with_restart_ig",
                       ]
    
    index4cutoff = ['distance_cutoff_uG']
    
    el.cut_records(normalize_indices=index4cutoff,smooth_length=300)  
    el.cut_and_smooth_normalize_records(normalize_indices=index4normalize,smooth_length=300)
    

    total_slice_num = 2 # edgelist 切成5段
    ordered_slice_num = 0 #选取要分析的一段, 可取 0,1,2,3,4,...total_slice_num-1
    index_min =0
    index_max =0
    length = len(el.records["degree"])
    
    
    if ordered_slice_num !=total_slice_num-1:
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
    
    m = int(np.power(len(index_data["if_rand"])/2,0.4)/2)
    index_data = mix2discrete_data(index_data,m)
        
    
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
    print(file_name+" 's final K is", progressive_removal(K))

    
    # sent_notifer()
    # useful_value,v,ave,(ci0,ci1),if_large_zero, multi =ep.casual_entropy('if_rand','distance',set(['indegree','degree']),index_data)
    # print(useful_value,v,ave,(ci0,ci1),if_large_zero,multi)
    pass

    
def analysis_index(el,total_slice_num,ordered_slice_num):
    index4normalize = ["indegree",
                       "degree",
                       "continuous_degree",
                       "k1sum_uG",
                       "k2sum",
                       "katz_uG",
                       "eigenvector_ig",
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
                       "add_time",
                       
                       #"resource_allocation_uG",
                       #"jaccard_coefficient_uG",
                       #"adamic_adar_index_uG",
                       "preferential_attachment_uG",
                       #"cn_soundarajan_hopcroft_uG",
                       #"ra_index_soundarajan_hopcroft_uG",
                       #"within_inter_cluster_uG",
                       "common_neighbor_centrality_uG",
                       "salton_uG",
                       "sorensen_uG",
                       "HPI_uG",
                       "HDI_uG",
                       "LHN1_uG",
                       "random_walk_with_restart_ig",
                       "efficiency_uG"
                       ]
    
    index4cutoff = ['shortest_path_length_uG','distance_cutoff_uG']
    index4addepsilon = ["resource_allocation_uG",
                    "jaccard_coefficient_uG","adamic_adar_index_uG","cn_soundarajan_hopcroft_uG",
                   "ra_index_soundarajan_hopcroft_uG","within_inter_cluster_uG"]
    
    smooth_length = 400
    el.cut_records(normalize_indices=index4cutoff,smooth_length=smooth_length)
    el.add_small_number_avoid_divide_by_zero(avoid_zero_indices=index4addepsilon,epsilon=0.01)    
    el.cut_and_smooth_normalize_records(normalize_indices=index4normalize+index4addepsilon,smooth_length=smooth_length)
    

    # total_slice_num = 2 # edgelist 切成5段
    # ordered_slice_num = 1 #选取要分析的一段, 可取 0,1,2,3,4,...total_slice_num-1
    index_min =0
    index_max =0
    length = len(el.records["degree"])
    
    
    if ordered_slice_num !=total_slice_num-1:
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
        
    m = int(np.power(len(index_data["if_rand"])/2,0.4)/2)
    index_data = mix2discrete_data(index_data,m)
    
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
    
    print(el.name +" 's final K is", progressive_removal(K))
    print(el.name +" 's final K is", progressive_removal(K),file=open("./experimental_results/output.txt", "a"))
    
    # sent_notifer()
    # useful_value,v,ave,(ci0,ci1),if_large_zero, multi =ep.casual_entropy('if_rand','distance',set(['indegree','degree']),index_data)
    # print(useful_value,v,ave,(ci0,ci1),if_large_zero,multi)
    pass



#file_name = './real_evolving_networks/hepph.txt'
#file_name = "./real_evolving_networks/hepth.txt"
#file_name = "./real_evolving_networks/sx-mathoverflow-a2q.txt"
#file_name = "./real_evolving_networks/sx-mathoverflow-c2a.txt"
#file_name = "./real_evolving_networks/sx-mathoverflow-c2q.txt"
#file_name = "./real_evolving_networks/wiki-talk-temporal.txt"
#file_name = "./real_evolving_networks/CollegeMsg.txt"
#file_name = "./real_evolving_networks/soc-sign-bitcoinalpha.txt"
# t0 = time.time()
el = network_evolve(file_name,max_edge_number=15000)
analysis_index(el,3,0)
analysis_index(el,3,1)
analysis_index(el,3,2)
# print("network evolve time is ",time.time()-t0)


