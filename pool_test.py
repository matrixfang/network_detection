import multiprocessing as mp
import threading as td
import time

def job(x):
    res = 0
    for i in range(x):
        res += i**2
    return res

def multicore_multiprocess():
    pool  = mp.Pool()
    for i in range(8):
        pool.apply_async(job,(100000000,))  #每次一个任务互相不影响
    pool.close()
    pool.join()

def multicore():
    pool  = mp.Pool(3)
    res = pool.map(job,range(10))#一次很多任务
    print(res)
    for i in range(8):
        job(100000000)  #每次一个任务互相不影响
    #print(res4.get())
    # multi_res = [pool.apply_async(job,(i,)) for i in range(10)]   
    # print([res.get() for res in multi_res])


if __name__ == "__main__":
    t0 = time.time()
    multicore()
    t1 = time.time()
    t2 = time.time()
    multicore_multiprocess()
    t3 = time.time()
    print("not multi finish time is", t1-t0)
    print("multi finish time is", t3-t2)