from evolving_mechanism import *

def network_evolve(file_name,max_edge_number=float('inf')):
    el = EdgeList() #build edgelist object
    index4normalize = ["indegree",
                       "degree",
                       "continuous_degree",
                       "k1sum_uG",
                       "k2sum",
                       "katz_uG",
                       #running time"eigenvector_uG",
                       "closeness_uG",
                       ## connected graph "information_uG",
                       ## running time "betweeness_centrality_uG",
                       ## connected graph "current_flow_betweenness_uG",
                       #"communicability_betweenness_uG",
                       "communicability_exp_uG",
                       
                       "harmonic_uG",
                       "pagerank_ig",
                       "h_index_uG",
                       "kshell_uG",
                       
                       
                       #"shortest_path_length_uG",
                       "resource_allocation_uG",
                       "jaccard_coefficient_uG",
                       "adamic_adar_index_uG",
                       "preferential_attachment_uG",
                       "cn_soundarajan_hopcroft_uG",
                       "ra_index_soundarajan_hopcroft_uG",
                       "within_inter_cluster_uG",
                       "common_neighbor_centrality_uG",
                       "salton_uG",
                       "sorensen_uG",
                       "HPI_uG",
                       "HDI_uG",
                       "LHN1_uG",
                       "average_commute_time_ig",
                       "cos_ig",
                       "random_walk_with_restart_ig",
                       "matrix_forest_index_ig",
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
    
    
     
    # el.cut_records(normalize_indices=index4cutoff,smooth_length=50)
    # el.cut_and_smooth_normalize_records(normalize_indices=index4normalize,smooth_length=50)

    sent_notifer()
    # el.load_records("BA_network_all.xlsx")

file_name = './real_evolving_networks/hepph.txt'
file_name = "./real_evolving_networks/hepth.txt"

network_evolve(file_name,max_edge_number=10000)