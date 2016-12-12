from __future__ import print_function
from pqt import PQTDecomposition
from helper_functions import *


def plg_alg(pd_edges, p_hat=0.01, pqt=None):
    """ Implementation of the PLG algorithm 
        Parameters:
        pd_edges - dictionary with pickup and delivery pairings
        p_hat - upperbound on the probability over a leaf of the pdf
        pqt - pass through a previously generated pdf

        Returns:
        dp_edges - lookup of computed delivery to pickup links
        """
    pickups = pd_edges.keys()
    deliveries = pd_edges.values()

    # If no pqt is passed then generate one from the given points
    if not pqt:
        pqt = PQTDecomposition().from_points(pickups, p_hat=p_hat)

    # Add all pickups to the tree
    pqt.add_points(pickups, 'p')

    dp_edges = {}
    
    cur_pickup = first_pickup = pickups.pop()
    if cur_pickup:
        # Find the leaf to remove the pickup
        pqt.enclosing_leaf(cur_pickup).content['p'].remove(cur_pickup)


    # While there are unvisited pickups
    while pickups:
        # Find the next delivery
        cur_delivery = pd_edges[cur_pickup]

        cur_leaf = pqt.enclosing_leaf(cur_delivery)

        # Use any local pickup if they exist
        if cur_leaf.content['p']:
            # Connect within leaf
            cur_pickup = cur_leaf.content['p'].pop()
            pickups.remove(cur_pickup)
        # Otherwise get a random unvisited pickup
        else:
            # Connect to any non-local pickup
            cur_pickup = pickups.pop()
            if cur_pickup:
                # Find the leaf to remove the pickup
                pqt.enclosing_leaf(cur_pickup).content['p'].remove(cur_pickup)

        # Add the edge
        dp_edges[cur_delivery] = cur_pickup

    # Add edge to make it a loop
    dp_edges[pd_edges[cur_pickup]] = first_pickup

    print(len(pqt.leaves))
    return dp_edges


def plg_test_1(n_pairs=50, verbose=False):
    pd_edges = gen_pd_edges(n_pairs=n_pairs)
    dp_edges = plg(pd_edges, p_hat=0.0025, pqt=None)
    
    if verbose:
        print_cycle(pd_edges, dp_edges)
 



if __name__ == "__main__":
    plg_test_1(n_pairs=10000, verbose=True)


    #decomp = PQTDecomposition().from_points(pts, p_hat=0.05, store=True)

    ##def pdf(x, y):
    ##    return 3 * (1 - x**2 - y**2)
    ##decomp = PQTDecomposition().from_pdf(pdf, p_hat=0.01, verbose=True)
    #print(decomp.enclosing_leaf([0.1,0.1]))
    #print(decomp)

   
