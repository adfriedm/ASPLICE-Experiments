from __future__ import print_function
import sys

from collections import defaultdict

from helper_functions import *
import numpy as np
import time
import cPickle as pickle

from splice import splice_alg


def save_res(res, file_pre="", dir_name=""):
    pickle.dump(res, open("{}{}_{}.pickle".format(dir_name,
        file_pre, res['experiment_time']), "wb"))

def splice_convergence_run(params=[(10,10),(20,10)]):

    res = defaultdict(list)
    res['experiment_time'] = time.strftime("%Y_%m_%d_%H_%M_%S", time.gmtime())

    # Counters
    sample = 0

    # Run experiment
    for (n_pairs, n_reps) in params:
        # 

        for rep in xrange(n_reps):
            print("{} pairs with rep ({}/{})".format(n_pairs, rep+1, n_reps))
            sys.stdout.flush()
            pd_edges = gen_pd_edges(n_pairs=n_pairs)
            sample += 1

            # Run experiment
            start_time = time.clock()
            dp_edges, cost = splice_alg(pd_edges)
            stop_time = time.clock()

            # Store the result
            res['sample'].append(sample)
            res['cost'].append(cost)
            res['time'].append(stop_time - start_time)
            res['n_pairs'].append(n_pairs)


    return res



def splice_convergence_exp1(res):
    import numpy as np
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.scatter(res['n_pairs'], res['cost'], alpha=0.5)
    
    ax.set_xlabel("num of pairs", fontsize=20)
    ax.set_ylabel("cost", fontsize=20)
    ax.set_title('SPLICE Cost vs Number of Pairs')
    
    ax.grid(True)
    fig.tight_layout()
    
    plt.show()

if __name__ == "__main__":

    res = splice_convergence_run([(int(n_samp), 2) for n_samp in
            np.logspace(1.0, 3.0, num=20)])
            
           
    save_res(res, file_pre="splice-conv",
            dir_name="results/splice-conv-test/")



#    map(lambda a,b,c,d: print("sample: {} cost: {:.3} time: {:.6} n_pairs {}".format(a,b,c,d)),
#        res['sample'],
#        res['cost'],
#        res['time'],
#        res['n_pairs'] )

    splice_convergence_exp1(res)

