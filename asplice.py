from __future__ import print_function
from pqt import PQTDecomposition
from helper_functions import *

import scipy.spatial.distance as dist

# PuLP Modeller for LP solver
import pulp


def asplice_alg(pd_edges, p_hat=0.01, pqt=None):
    pickups = pd_edges.keys()
    deliveries = pd_edges.values()
    # If no pqt is passed then generate
    if not pqt:
        pqt = PQTDecomposition().from_points(pickups, p_hat=p_hat)

    # Add all pickups and deliveries to the tree
    pqt.add_points(pickups, 'p')
    pqt.add_points(deliveries, 'd')

    # Create a d->p edge dictionary
    dp_edges = {}

    # For storing leaves with excess after locally connecting
    excess_p_leaves = []
    excess_d_leaves = []

    # For every leaf in the decomposition
    for leaf in pqt.leaves:
        # If there are excess pickups and deliveries then connect them
        while leaf.content['p'] and leaf.content['d']:
            pickup = leaf.content['p'].pop()
            delivery = leaf.content['d'].pop()
            dp_edges[delivery] = pickup
        else:
            # Store a leaf if it has excess
            if len(leaf.content['p']) > 0:
                excess_p_leaves.append(leaf)
            elif len(leaf.content['d']) > 0:
                excess_d_leaves.append(leaf)

    # Note: We are only left with excess of either pickups or deliveries,
    # the total of which can be shown to be O(n^.5).

    # Now form the transportation problem on the leaves' excess
    # Compute distance matrix
    cost_mat = dist.cdist(
            map(lambda leaf: leaf.center(), excess_d_leaves),
            map(lambda leaf: leaf.center(), excess_p_leaves))
    # Map edges to cost
    cost_dict = {(leaf_d, leaf_p): cost_mat[i][j]
            for i, leaf_d in enumerate(excess_d_leaves)
            for j, leaf_p in enumerate(excess_p_leaves)}
    # Create all possible excess edges
    dp_intra = [(leaf_d, leaf_p) for leaf_d in excess_d_leaves
            for leaf_p in excess_p_leaves]

    # Setup LP model
    model = pulp.LpProblem("LP Problem", pulp.LpMinimize)
    # Pulp variables
    x = pulp.LpVariable.dicts("edgeweight", (excess_d_leaves,
        excess_p_leaves), 0, None, pulp.LpInteger)
    # Objective function
    model += pulp.lpSum([x[leaf_d][leaf_p]*cost_dict[(leaf_d,leaf_p)]
        for (leaf_d,leaf_p) in dp_intra])
    # Constraints
    for leaf_d in excess_d_leaves:
        model += pulp.lpSum([x[leaf_d][leaf_p] for leaf_p in
            excess_p_leaves]) == len(leaf_d.content['d']), \
            "contrain delivery excess for {}".format(leaf_d)
    for leaf_p in excess_p_leaves:
        model += pulp.lpSum([x[leaf_d][leaf_p] for leaf_d in
            excess_d_leaves]) == len(leaf_p.content['p']), \
            "contrain pickup excess for {}".format(leaf_p)

    # Solve the LP
    status = model.solve()

    # Connect greedily from leaves with excess delivery to those with
    # excess pickups
    for dp in dp_intra:
        leaf_d,leaf_p = dp
        for i in xrange(int(x[leaf_d][leaf_p].value())):
            dp_edges[leaf_d.content['d'].pop()] = leaf_p.content['p'].pop()
    cost = merge_tours(pd_edges, dp_edges, dist_f=dist.euclidean)

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


def asplice_test_1(n_pairs=100, verbose=False):
    pd_edges = gen_pd_edges(n_pairs=n_pairs)
    dp_edges, cost = asplice_alg(pd_edges, p_hat=0.0025, pqt=None)
    
    print("Cost: {}".format(cost))

    if verbose:
        print_cycle(pd_edges, dp_edges)


if __name__ == '__main__':
    asplice_test_1(n_pairs=1000, verbose=False)

