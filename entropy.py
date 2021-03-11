from npeet import entropy_estimators as ee
import numpy as np
from collections import Counter
from scipy import stats


def discrete_entropy(dxs):
    dic = dict(Counter(dxs))
    sum = np.sum(list(dic.values()))
    ent = 0.0
    for k in dic:
        p = float(dic[k])/sum
        ent += - p * np.log2(p)
    return ent

def condition_entropy(dxs,dys):
    H = 0.0
    fvalues2 = {}
    values1 = {}
    for ind in range(len(dys)):
        iterm = dys[ind]
        if iterm in values1:
            values1[iterm].append(dxs[ind])
            fvalues2[iterm]+=1
        else:
            values1[iterm]=[]
            values1[iterm].append(dxs[ind])
            fvalues2[iterm] = 1
    for key in fvalues2:
        p = float(fvalues2[key])/len(dys)
        H += p * discrete_entropy(values1[key])
    return H

def slice(xs,m):
    #devide into m spaces
    ls = sorted(xs)
    npers = int(len(xs)/m)
    dx = []
    for k in range(0, m):
        dx.append(float(ls[k * npers]))
    dx.append(float(ls[len(ls)-1]+1.0))
    # print(dx)
    # print(len(dx))

    def rep(v):
        for i in range(0, len(dx)-1):
            min = dx[i]
            max = dx[i + 1]
            if v>=min and v< max: #不确定是slice是否对
                return i
        return len(dx)-1
        
    return np.fromiter(map(rep, xs),dtype = int)

def slice_by_random(xs,x_random,m):
    #devide into m spaces
    #按照x_random 的分布切分 xs 使其变成离散值
    ls = sorted(x_random)
    npers = int(len(x_random)/m)
    dx = []
    for k in range(0, m):
        dx.append(float(ls[k * npers]))
    dx.append(float(ls[len(ls)-1]+1.0))
    # print(dx)
    # print(len(dx))

    def rep(v):
        for i in range(0, len(dx)-1):
            min = dx[i]
            max = dx[i + 1]
            if v>=min and v< max: #不确定是slice是否对
                return i
        return len(dx)-1
        
    return np.fromiter(map(rep, xs),dtype = int)



def cmiddcut(x,y,z,slice_num=20):
    #under condition z
    x_vec = np.transpose(np.array([x]))
    y_vec = np.transpose(np.array([y]))
    z_vec = np.transpose(np.array([slice(z, slice_num)]))
    return ee.cmidd(x_vec,y_vec,z_vec)

def cmidcutd(x,y,z,slice_num=20):
    #under condition z
    x_vec = np.transpose(np.array([x]))
    y_vec = np.transpose(np.array([slice(y,slice_num)]))
    z_vec = np.transpose(np.array([z]))
    return ee.cmidd(x_vec,y_vec,z_vec)

def midcut(x,y,slice_num=20,base=2):
    x_vec = np.transpose(np.array([x]))
    y_vec = np.transpose(np.array([slice(y, slice_num)]))
    return ee.midd(x_vec,y_vec)
def midd(x,y,slice_num=20,base=2):
    x_vec = np.transpose(np.array([x]))
    y_vec = np.transpose(np.array([y]))
    return ee.midd(x_vec,y_vec)
def shuffle_test_dcut(x,y,slice_num=20,base=2):
    x_vec = np.transpose(np.array([x]))
    y_vec = np.transpose(np.array([slice(y, slice_num)]))
    # z_vec = np.transpose(z)
    print(ee.shuffle_test(ee.midd,x_vec,y_vec,base=base))


def midc(x,y,k=10,base=2):
    x_vec = np.array([[s] for s in x ])
    y_vec = np.array([[s] for s in y ])
    # print(x_vec)
    # print(y_vec)
    return ee.midc(x_vec,y_vec,k=k,base=base)

def cmiddc(x,y,z,k =10,base=2):
    #entropyd(xz, base) + entropyd(yz, base) - entropyd(xyz, base) - entropyd(z, base)
    x_vec = np.array([[s] for s in x ]) # d
    y_vec = np.array([[s] for s in y ]) #d
    z_vec = np.array([[s] for s in z ]) # c
    xy_vec = np.c_[x_vec, y_vec]
     
    cmiddc = ee.midd(x_vec,y_vec,base=base) - (ee.midc(y_vec,z_vec,k=k,base=base,warning=False)+ee.midc(x_vec,z_vec,k=k,base=base,warning=False)-ee.midc(xy_vec,z_vec,k=k,base=base,warning=False))
    cmidcd = ee.midc(x_vec,z_vec,k=k,base=base,warning=False) - (ee.midc(y_vec,z_vec,k=k,base=base,warning=False)+ee.midc(x_vec,z_vec,k=k,base=base,warning=False)-ee.midc(xy_vec,z_vec,k=k,base=base,warning=False))
    
    
    return cmiddc,cmidcd

def casual_entropy(i,j,K,data):    

    x = data[i]
    y = data[j]
    if len(K) ==0:
        return casual_entropy_empty(i,j,K,data)
    
    #slice_num = int(np.power(len(x),1.0/2)/2)
    slice_num = int(np.power(len(x),0.4)/2)
      
    if x.dtype =='float64':
        x = slice(x,slice_num)
        x_vec = np.array([[s] for s in x ]) 
    else:
        x_vec = np.array([[s] for s in x ]) 
        
    if y.dtype =="float64":
        y = slice(y,slice_num)
        y_vec = np.array([[s] for s in y ]) 
    else:
        y_vec = np.array([[s] for s in y ]) 
        print("index j "+j +"is discrete")
    
    
    z_all = []    
    for k in K:
        z = data[k]
        if z.dtype =="float64":
            z = slice(z,slice_num)
            z_vec = np.array([[s] for s in z ]) 
        else:
            z_vec = np.array([[s] for s in z ]) 
        z_all.append(z_vec)
    z_combine = np.c_[tuple(z_all)]


    
    #x_clone = np.copy(x_vec)
    y_clone = np.copy(y_vec)
    ns = 200
    ci = 0.95
    outputs = []
    #outputs2 = []
    for i in range(ns):
        np.random.shuffle(y_clone)
        outputs.append(ee.cmidd(x_vec,y_clone,z_combine,base=2))
        # outputs2.append(ee.midd(x_clone,y_vec,base=2))
    outputs.sort()
    # outputs2.sort()

    v = ee.cmidd(x_vec,y_vec,z_combine,base=2)
    ave = np.mean(outputs)
    ci0 = outputs[int((1. - ci) / 2 * ns)] 
    ci1 = outputs[int((1. + ci) / 2 * ns)] 
     
    if v > ci1:
        if_large_zero = True
    else:
        if_large_zero = False
        
    
    n =200
    useful_result = if_large_zero * v
    std_modified = np.sqrt((n-1)/n)*np.std(outputs)
    # multi = abs(v-ave)/np.std(outputs)
    multi = 0
    
    #(statistic, pvalue) = stats.ttest_ind_from_stats(mean1=ave, std1=std_modified, nobs1=200, mean2=v, std2=0, nobs2=2,equal_val=False)
    # res = stats.ttest_1samp(np.array(outputs),[ave,0])
    statistic = 1
    
    
    print("index "+j+" function: casual_entropy, the length of data",len(x),"the slice number is", slice_num, "useful value is ",useful_result,"multi sigma is ", multi)
   # print("statistic, pvalue",statistic,pvalue)  
    
    return useful_result,v,ave,(ci0,ci1),if_large_zero,abs(statistic)*if_large_zero
    
def casual_entropy_empty(i,j,K,data):    

    x = data[i]
    y = data[j]
    
    slice_num = int(np.power(len(x),0.4)/2)
 
    if x.dtype =='float64':
        x = slice(x,slice_num)
        x_vec = np.array([[s] for s in x ]) 
    else:
        x_vec = np.array([[s] for s in x ]) 
        
    if y.dtype =="float64":
        y = slice(y,slice_num)
        y_vec = np.array([[s] for s in y ]) 
    else:
        y_vec = np.array([[s] for s in y ]) 
    
        

    # x_clone = np.copy(x_vec)
    y_clone = np.copy(y_vec)
    ns = 200
    ci = 0.95
    outputs = []
    outputs2 = []
    for i in range(ns):
        np.random.shuffle(y_clone)
        outputs.append(ee.midd(x_vec,y_clone,base=2))
        # outputs2.append(ee.midd(x_clone,y_vec,base=2))
    outputs.sort()
    # outputs2.sort()

    v = ee.midd(x_vec,y_vec,base=2)
    ave = np.mean(outputs)
    ci0 = outputs[int((1. - ci) / 2 * ns)] 
    ci1 = outputs[int((1. + ci) / 2 * ns)] 
     
    if v > ci1:
        if_large_zero = True
    else:
        if_large_zero = False
    
    useful_result = if_large_zero * v
    # multi = abs(v-ave)/np.std(outputs)
    multi = 0
    
    n=200  #check if statistical significant
    # std_modified = np.sqrt((n-1)/n)*np.std(outputs)
    # (statistic, pvalue) = stats.ttest_ind_from_stats(mean1=ave, std1=std_modified, nobs1=200, mean2=v, std2=0, nobs2=1)
    statistic =1
    
    print("index "+j+" function: casual_entropy, the length of data",len(x),"the slice number is", slice_num, "useful value is ",useful_result,"multi sigma is ", multi) 
    # print("statistic, pvalue",statistic,pvalue)  
    
    return useful_result, v,ave,(ci0,ci1),if_large_zero,abs(statistic)*if_large_zero
 


def shuffle_midd(x,y,k =10,base=2):
    print(ee.shuffle_test(ee.midd,x,y,base=base,k = k))
    pass


def entropy_package(x,y,z):
    # x_d , y_d, z_c
    
    a = cmiddcut(x,y,z)
    b = cmidcutd(x,z,y)
    c = midcut(x, z,slice_num=20)
    d = midd(x,y)
    return a,b,c,d
    
    
    pass
def cmiddd():
    pass

def cmidcd():
    pass

if __name__ == '__main__':
    x = np.array([1,1,1,1,2,2,2,2])
    y = np.array([1.1,1.01,2.1,2-0.1,1+0.02,1-0.05,2+0.1,2-0.2])
    z = np.array([1,1,1,2,2,2,2,2])
    print(x.reshape((8,1)))
    
    print(ee.entropy(x.reshape(len(x),1),base=2))
    print(slice(y,2))
    print(ee.midd(x,y,base=2))
    print(ee.shuffle_test(ee.midd,x,y,base=2))


          
          
          