from __future__ import print_function
from collections import defaultdict

from helper_functions import *
import numpy as np
import time
import cPickle as pickle

from splice import splice_alg


def save_res(res, file_pre="", dir_name=""):
    pickle.dump(res, open("{}{}_{}.pickle".format(dir_name,
        file_pre, res['experiment_time']), "wb"))

def splice_convergence_test(n_reps=20, n_pairs_l=20, n_pairs_h=23):
    res = defaultdict(list)
    res['experiment_time'] = time.strftime("%Y_%m_%d_%H_%M_%S", time.gmtime())

    # Counters
    sample = 0

    # Run experiment
    for n_pairs in xrange(20, 23):
        for rep in xrange(20):
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




if __name__ == "__main__":

    res = splice_convergence_test(n_reps=20, n_pairs_l=20, n_pairs_h=23)

    save_res(res, file_pre="splice-conv",
            dir_name="results/splice-conv-test/")

    map(lambda a,b,c,d: print("sample: {} cost: {:.3} time: {:.6} n_pairs {}".format(a,b,c,d)),
        res['sample'],
        res['cost'],
        res['time'],
        res['n_pairs'] )



