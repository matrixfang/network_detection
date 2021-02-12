import networkx as nx
from collections import Counter
import evolving_mechanism as em
import matplotlib.pyplot as plt
import numpy as np

class pafit(object):
    def __init__(self):
        super().__init__()
        self.edge_list = []
        self.m = []
        self.n = []
        self.g = None
    def load_edge_list(self,el):
        self.edge_list = el
    def evolve(self):
        g = nx.DiGraph()
        for t in range(100):
            edge = self.edge_list[t]
            g.add_edge(edge[0],edge[1])
        for t in range(100,len(self.edge_list)):
            edge = self.edge_list[t]
            ni = edge[0]
            nj = edge[1]
            if nj in g.nodes():
                self.m.append({g.degree(nj):1})
                self.n.append(dict(Counter(dict(g.degree()).values())))
                g.add_edge(ni,nj)
            else:
                self.m.append({0:1})
                self.n.append(dict(Counter(dict(g.degree()).values())))
                g.add_edge(ni,nj)
                
        self.g = g
    
    def pafit(self):
        K = max(dict(self.g.degree()).values())
        a = []
        for k in range(K+1):
            a.append(1.0)
        b = np.array(a)
        
        def l(a):
            K = len(a)
            s = 0
            for t in range(len(self.m)):
                s2 = 0
                s1 =0
                for k in range(K+1):
                    if (a[k]) == 0:
                        s1 += 0
                    else:
                        s1 += self.m[t][k] * np.log(a[k])
                    s2 += self.n[t][k] * a[k]
                s+=s1- sum((self.m[t]).values()) *np.log(s2)
                 
            return s
        
        def iteration(a):
            over_k = {}
            under_k = {}
            na_t = []
            
            for t in range(len(self.m)):
                m = self.m[t]
                for key in m:
                    if key in over_k:
                        over_k[key] += m[key]
                    else:
                        over_k[key] = 0
            # print("over_k",over_k)
                    
            for t in range(len(self.m)):
                na_t.append(0)
                for k in self.n[t]:
                    na_t[t] += a[k] * self.n[t][k]
            # print("na_t",na_t)
            
            for t in range(len(self.m)):
                for k in self.n[t]:
                    if k in under_k:
                        under_k[k] += self.n[t][k]* sum((self.m[t]).values())/na_t[t]
                    else:
                        under_k[k] = 0
            
            # print("under_k",under_k)
            
            for k in range(K+1):
                if k in over_k and under_k[k] != 0:
                    a[k] = over_k[k]/under_k[k]
                else:
                    a[k] = 0
            # print("a",a)
            scale = a[1]
            a /= a.sum()
            # for k in range(K+1):
            #     a[k] = a[k]/scale
        
        for i in range(20):
            iteration(b)
            # print(l(a))
        return a
            
def test():
    p = pafit()
    el = em.EdgeList()
    g = el.read_BA(0.8,1000)
    # g = el.read_HK(1.2,0.0,1000)
    # edge_list = [(1,2),(3,2),(4,1),(5,4),(6,7)]
    p.load_edge_list(el.edge_list)
    p.evolve()
    a = p.pafit()
    plt.plot(a)
    plt.show()
    
    pass
if __name__ == "__main__":
    test()
    