from __future__ import print_function
#import numpy as np
#from pqt import PQTDecomposition, PQTNode

from random import random

def print_cycle(pd_edges, dp_edges):
    # Print the cycle
    first_pickup = cur_pickup = next(key for key in pd_edges)
    print("   P ({:.3},{:.3}) ".format(*cur_pickup))

    # Do while construction
    while True:
        cur_delivery = pd_edges[cur_pickup]
        cur_pickup = dp_edges[cur_delivery]
        print("-> D ({:.4},{:.4}) ".format(*cur_delivery))
        # Break out if back at beginning
        if cur_pickup == first_pickup:
            print("-> Loop")
            break
        else:
            print("-> P ({:.4},{:.4})".format(*cur_pickup))


def gen_pd_edges(n_pairs=50):
    return {(random(),random()): (random(),random()) \
                for i in xrange(n_pairs)}
  
