import numpy as np
import matplotlib.pyplot as plt
import time
from evolving_mechanism import *
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


def test_networkx():
    el = EdgeList()
    g = el.read_BA(1.0)
    print(g.nodes,type(g.nodes()))


test_container()