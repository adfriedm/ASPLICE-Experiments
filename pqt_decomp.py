from __future__ import print_function
import numpy as np

class PQTNode:
    """PQT Node class"""
    def __init__(self, coords=[[0., 1.], [0., 1.]]):
        self.children = []
        self.coords = coords
        self.contents = []
        self.p = 0.

    def __str__(self):
        return "Node {}x{} \t{} points {} children".format(self.coords[0], self.coords[1], len(self.contents), len(self.children))

    def __repr__(self):
        return "PQTNode({}, {})".format(self.coords[0], self.coords[1])

    def split(self):
        x0, x1 = self.coords[0]
        y0, y1 = self.coords[1]

        xc, yc = 0.5*(x0+x1), 0.5*(y0+y1)

        # Add subcoordinates
        self.children = [
                PQTNode([[x0,xc],[y0,yc]]),
                PQTNode([[xc,x1],[y0,yc]]),
                PQTNode([[xc,x1],[yc,y1]]),
                PQTNode([[x0,xc],[yc,y1]])
                ]
        return self.children

    def contains(self, coord=[0.5,0.5]):
        x0, x1 = self.coords[0]
        y0, y1 = self.coords[1]

        return x0 <= coord[0] < x1 \
           and y0 <= coord[1] < y1


class PQTDecomposition:
    """PQT Decomposition data structure class"""
    def __init__(self):
        self.root = PQTNode()
        self.leaves = []

    def from_points(self, points=[], p_hat=0.1, store_pts=True):
        n_pts = float(len(points))
        # Check that atoms do not have probability higher than p_hat, if they
        # are then we set p_hat to something slightly bigger than this value.
        atom_p = 1./n_pts
        self.p_hat = atom_p if (atom_p > p_hat) else p_hat


        def gen_pqt(node, pts):
            node.p = len(pts)/n_pts

            # The first condition is the subpartitioning rule for a pqt.
            if node.p >= p_hat and len(pts) > 1:
                # Add children to the current node
                node.split()
                # For each new node, generate from all points that fall inside
                # the cell
                for child in node.children:
                    gen_pqt(child, [pt for pt in pts if child.contains(pt)])
                return

            # Otherwise the node is a leaf, so add it
            self.leaves.append(node)
            # Store points if set
            if store_pts:
                node.contents = pts

        # Start recursion through the root node
        gen_pqt(self.root, points)
        return self

    def __ref__(self):
        return "PQTDecomposition()"

    def __str__(self):
        print_str = ""

        # Store node, depth data on stack. Work through tree depth first
        node_stack = [(self.root, 0)]
        # If there are things on the stack
        while node_stack:
            node, depth = node_stack.pop()
            i = None
            for i in xrange(depth):
                print_str += "  "
            else:
                if i is not None:
                    print_str += "- "

            print_str += str(node) + "\n"

            # If the node has children then process them next on the stack
            for child in node.children:
                node_stack.append((child,depth+1))
        return print_str

        

if __name__ == "__main__":
    import random
    n_pts = 1000 
    pts = [(random.random(), random.random()) for i in xrange(n_pts)]
    decomp = PQTDecomposition().from_points(pts, p_hat=0.05)
    
    map(print, decomp.leaves)

