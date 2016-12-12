from __future__ import print_function
from pqt import PQTDecomposition
from helper_functions import *

def asplice_alg(pd_edges, p_hat=0.01, pqt=None):
    pickups = pd_edges.keys()
    deliveries = pd_edges.values()

    # If no pqt is passed then generate
    if not pqt:
        pqt = PQTDecomposition().from_points(pickups, p_hat=p_hat)

    # Add all pickups and deliveries to the tree
    pqt.add_points(pickups, 'p')
    pqt.add_points(deliveries, 'd')


    dp_edges = {}

    # For every leaf in the decomposition
    for leaf in pqt.leaves:
        # If there are excess pickups and deliveries then connect them
        while leaf.content['p'] and leaf.content['d']:
            pickup = leaf.content['p'].pop()
            delivery = leaf.content['d'].pop()
            dp_edges[delivery] = pickup
    
    # We are only left with excess of either pickups or deliveries,
    # the total of which can be shown to be O(n^.5).
    # Now form the transportation problem on the leaves of the
    # decomposition


    # Connect greedily from

