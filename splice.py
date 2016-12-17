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

def splice_alg(pd_edges):
    """ Implementation of the SPLICE algorithm
    """
    pickups = pd_edges.keys()
    deliveries = pd_edges.values()

    # Compute Hungarian matching indices
    c_mat = dist.cdist(deliveries, pickups)
    # Connect dp_edges according to bipartite matching
    dp_edges = {deliveries[d_ind]: pickups[p_ind] \
                for (d_ind, p_ind) in Munkres().compute(c_mat)}
    # Merge generated tours
    cost = merge_tours(pd_edges, dp_edges)

    return dp_edges, cost


def merge_tours(pd_edges, dp_edges, dist_f=dist.euclidean):
    cost = 0.
    pickups = pd_edges.keys()
    deliveries = pd_edges.values()

    # Setup beginning of merge
    cur_p = pickups.pop()
    tour_p = cur_p
    start_p = cur_p

    # While there are remaining pickups
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

    return cost
    




if __name__ == "__main__":
    pass
