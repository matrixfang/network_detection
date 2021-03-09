from evolving_mechanism import *
import os
el = EdgeList()

file_name = "./real_evolving_networks/soc-sign-bitcoinalpha.csv"
el.read_edgelist_file(file_name)

with open("./soc.txt","w") as f:
    for e in el.edge_list:
        f.write(str(e[0])+" "+str(e[1])+" "+str(e[2])+os.linesep)