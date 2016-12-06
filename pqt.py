from __future__ import print_function
import numpy as np


class PQTNode:
    """PQT Node class"""
    def __init__(self, bounds=[[0., 1.], [0., 1.]]):
        self.children = []
        self.bounds = bounds 
        self.content = []
        self.p = 0.

    def __str__(self):
        return "[{:.3},{:.3}]x[{:.3},{:.3}] ".format(self.bounds[0][0],
                                                     self.bounds[0][1],
                                                     self.bounds[1][0],
                                                     self.bounds[1][1]) \
              + "{} pts {} chldrn {:.3} prb".format(len(self.content), len(self.children), self.p)

    def __repr__(self):
        return "PQTNode({}, {})".format(self.bounds[0], self.bounds[1])

    def split(self):
        x0, x1 = self.bounds[0]
        y0, y1 = self.bounds[1]

        xc, yc = 0.5*(x0+x1), 0.5*(y0+y1)

        # Add subcoordinates
        self.children = [
                PQTNode([[x0,xc],[y0,yc]]),
                PQTNode([[xc,x1],[y0,yc]]),
                PQTNode([[xc,x1],[yc,y1]]),
                PQTNode([[x0,xc],[yc,y1]])
                ]
        return self.children

    def encloses(self, coord=[0.5,0.5]):
        x0, x1 = self.bounds[0]
        y0, y1 = self.bounds[1]

        return x0 <= coord[0] < x1 \
           and y0 <= coord[1] < y1

    def contains(self, coord=[0.1, 0.1]):
        return coord in self.content


class PQTDecomposition:
    """PQT Decomposition data structure class"""
    def __init__(self):
        self.root = PQTNode()
        self.leaves = []

    def from_points(self, points=[], p_hat=0.1, store=False):
        n_pts = float(len(points))
        # Check that atoms do not have probability higher than p_hat, if they
        # are then we set p_hat to the probability of an atom.
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
                    gen_pqt(child, [pt for pt in pts if child.encloses(pt)])
            else:
                # Otherwise the node is a leaf, so add it
                self.leaves.append(node)
                # Store points if set
                if store:
                    node.content = pts

        # Start recursion through the root node
        gen_pqt(self.root, points)
        return self

    def from_pdf(self, pdf, p_hat):
        from scipy.integrate import nquad

        self.p_hat = p_hat

        def gen_pqt(node):
            # Compute the probability over the cell
            node.p,_ = nquad(pdf, node.bounds)

            # If the probability is too high then split the cell and generate
            # sub-trees
            if node.p >= p_hat:
                node.split()
                for child in node.children:
                    gen_pqt(child)
            else:
                # Otherwise the node is a leaf
                self.leaves.append(node)

        gen_pqt(self.root)
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

    def enclosing_leaf(self, coords):
        def _get_leaf(node):
            # Check all children (if any)
            for child in node.children:
                # Search down branch if contains coord
                if child.encloses(coords):
                    return _get_leaf(child)
            return node
        # Check if the point is enclosed by the pqt
        if self.root.encloses(coords):
            return _get_leaf(self.root)
        return None

    def add_point(self, coord):
        leaf = self.enclosing_leaf(coord)
        if not leaf:
            return False
        leaf.content.append(coord)
        return True
        
    def add_points(self, coords):
        # Cant be run in parallel because node tree is not
        # serializable
        return map(self.add_point, coords)


if __name__ == "__main__":
    pass
    #decomp = PQTDecomposition().from_points(pts, p_hat=0.05, store=True)
    
    ##def pdf(x, y):
    ##    return 3 * (1 - x**2 - y**2)
    ##decomp = PQTDecomposition().from_pdf(pdf, p_hat=0.01, verbose=True)
    #print(decomp.enclosing_leaf([0.1,0.1]))
    #print(decomp)