#!/usr/bin/env sage

'''
USAGE Ensure SAGE_ROOT is stored in your PATH, and run

./skein-dimensions.sage n

where n is an integer giving the number of levels to search through in a
triangular shell (see description), to estimate the dimension of the skein
module of T^3.

DESCRIPTION A script to estimate the dimension of the skein module of the
3-torus. Works by increasing the size of a triangular lattice shell: all points
in a box above the line y = -x. E.g:

                        x x x x x
                          x x x x
                            0 x x
                              x x
                                x

As shell_level -> infty this is a generating set. The program finds linear
relations between lattice points at each level, and then does row reduction to
compute the dimension. If the dimensions stabilize as n grows, we expect this
gives an upper bound on the dimension of the skein module.

Where lattice point (a, b) corresponds to generator m_(a, b), the linear
relations are given by

(q^{-bc} - q^{-ad})m_{a+c, b+d} + (q^{bc} - q^{ad})m_{c-d, b-d} = 0

for lattice points (a, b) and (c, d).

TODO: Implement a version for tori twisted by an element in SL_2(Z).
'''

import sys
from sage.all import *

def order_lrtb(shell_level):
    '''
    Returns a dictionary giving an order to the lattice points in a triangular
    shell.
    Ordering given by: (a, b) < (c, d) if a < c; (a, b) < (a, c) if b > c.
    That is, ordering goes left-right, top-bottom (lrtb) through the triangular
    shell.
    '''
    order_dict = {}
    place = 0 # Place of the current lattice point in the order (to increment)
    # Incrementally loop over a, then decrement through b until we reach y = -x.
    for a in range(-1*shell_level, shell_level+1):
        b = shell_level
        while (a + b) >= 0:
            # Populate dict entry, increment place for next point,
            # decrement y coord.
            order_dict.update({(a, b): place})
            place += 1
            b -= 1
    return order_dict

def get_relations(shell_level, order_func):
    '''
    Returns a list of linear relations between lattice points for a specified
    shell level. Each relation is a list representing a relation vector.

    Requires an integer shell_level, and a function order_func : int -> dict
    which should produce a dictionary with keys being lattice points in a
    certain shell level, and values being their position in some sequential
    ordering specified by the function (this is required to store the relations
    as vectors).

    Performs a double loop through the lattice, obtains the relation for each
    pair of points, and discards trivial or out-of-range relations.
    '''

    relations = []
    N = (2*shell_level + 1)*(shell_level + 1) # Total number of lattice points.
    ordering = order_func(shell_level) # Dictionary giving the points an order.

    # Incrementally loop over a, and decrement through b until we reach y = -x.
    # This loops through all lattice points (a, b) in this shell level.
    for a in range(-1*shell_level, shell_level+1):
        b = shell_level
        while (a + b) >= 0:
            # Similarly, loop over all points (c, d) in the shell at this level.
            for c in range(-1*shell_level, shell_level + 1):
                d = shell_level
                while (c + d) >= 0:
                    # Rewrite indices for the generators appearing in the linear
                    # relation coming from points ((a, b). (c, d)).
                    r = a + c
                    s = b + d
                    t = a - c
                    u = b - d
                    # Check the relations are not out of range.
                    if r + s >= 0 and t + u >= 0 and r <= shell_level and s <= shell_level and t <= shell_level and u <= shell_level:
                        # Compute the coefficients in the relations.
                        Q_1 = q**(-b*c) - q**(-a*d)
                        Q_2 = q**(b*c) - q**(a*d)
                        # Check the relations are not trivial.
                        if Q_1 != 0 and Q_2 != 0:
                            # Update the list of relations with a vector.
                            x = ordering[(r, s)]
                            y = ordering[(t, u)]
                            rel = [0]*N
                            rel[x] = Q_1
                            rel[y] = Q_2
                            relations.append(rel)
                    d -= 1
            b -= 1
    return relations

# Handle incorrect number of arguments.
if len(sys.argv) != 2:
    print("Give an integer n, to compute dimension estimates for n levels of triangular shells.")
    sys.exit(1)

n = sage_eval(sys.argv[1]) # Number of levels of shell to compute for.
q = var('q') # Declare an indeterminate q.

for shell_level in range(n+1):
    N = (2*shell_level + 1)*(shell_level + 1) # For each shell level, compute #{lattice points}.
    print("Calculating relations for level %d (%d lattice points) ..." % (shell_level, N))
    relations = get_relations(shell_level, order_lrtb)
    print("Found %d (non-independent) relations. Reducing ..." % len(relations))
    # Form a relation matrix, compute its rank; the dimension estimate is the
    # co-rank.
    A = matrix(CC['q'].fraction_field(), relations)
    relations_rank = A.rank()
    dim_estimate = N - relations_rank
    print("Dimension estimate for level %d: %d.\n" % (shell_level, dim_estimate))
