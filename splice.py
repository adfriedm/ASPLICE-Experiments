from __future__ import print_function
from helper_functions import *
from munkres import Munkres
import numpy as np
import scipy.spatial.distance as dist

import time

import cPickle as pickle


from collections import defaultdict

# Returns formatted string from list of tuples
def pt_str(pt):
    return "({0:.4},{1:.4})".format(*pt)
def pts_str(pts):
    return ", ".join(map(lambda pt: pt_str(pt), pts))

def splice_alg(pd_edges, p_hat=0.01):
    """ Implementation of the SPLICE algorithm
    """
    start_time = time.clock()
    pickups = pd_edges.keys()
    deliveries = pd_edges.values()

    cost = 0.

    # Compute Hungarian matching indices
    c_mat = dist.cdist(deliveries, pickups)
    # Connect dp_edges according to bipartite matching
    dp_edges = {deliveries[d_ind]: pickups[p_ind] \
                for (d_ind, p_ind) in Munkres().compute(c_mat)}

    cur_p = pickups.pop()
    tour_p = cur_p
    start_p = cur_p

    #import pdb; pdb.set_trace()
    while pickups:
        # Follow through the d->p chain
        cur_d = pd_edges[cur_p]
        next_p = dp_edges[cur_d]

        # If tour finished
        if next_p == tour_p:
            # Chose any random unvisited pickup
            next_p = pickups.pop()
            # Start new tour
            tour_p = next_p
            # Change dp to connect to new pickups
            dp_edges[cur_d] = next_p
        else:
            # Mark pickup as visited
            pickups.remove(next_p)
        cur_p = next_p

    dp_edges[pd_edges[cur_p]] = start_p

    # Sum over all pd and dp edge costs
    cost += reduce(lambda a, b: a + dist.euclidean(b,pd_edges[b]),
            pd_edges, 0)
    cost += reduce(lambda a, b: a + dist.euclidean(b,dp_edges[b]),
            dp_edges, 0)

    stop_time = time.clock()

    return dp_edges, cost, (stop_time - start_time)
    
    



def splice_convergence_test():
    res = defaultdict(list)
    res['exp_time'] = time.strftime("%Y_%m_%D_%H_%M_%S", time.gmtime())
    # Parameters
    n_pairs_l = 20
    n_pairs_h = 23
    n_reps = 20

    # Counters
    sample = 0

    # Run experiment
    for n_pairs in xrange(20, 23):
        for rep in xrange(20):
            pd_edges = gen_pd_edges(n_pairs=n_pairs)
            sample += 1
            dp_edges, cost, _time = splice_alg(pd_edges)
            res['sample'].append(sample)
            res['cost'].append(cost)
            res['time'].append(_time)
            res['n_pairs'].append(n_pairs)

    map(lambda a,b,c,d: print("sample: {} cost: {:.3} time: {:.6} n_pairs {}".format(a,b,c,d)),
        res['sample'],
        res['cost'],
        res['time'],
        res['n_pairs'] )

    # Pickle the results and store them 
    sdir_name = "splice-conv-test/"
    sfile_name = "splice-res_{}.pickle".format(res['exp_time'])

    pickle.dump(res, open(sdir_name+sfile_name, "wb"))


if __name__ == "__main__":

    splice_convergence_test()

