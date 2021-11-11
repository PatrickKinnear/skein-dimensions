#!/usr/bin/env sage

'''
USAGE Ensure SAGE_ROOT is stored in your PATH, and run

./skein-dimensions.sage

You will be prompted for an SL_2(Z)-matrix and an integer n, where n givesthe
number of levels to search through in a triangular shell (see description), to
estimate the dimension of the skein module of T^2 x S^1 twisted by the specified
matrix.

DESCRIPTION A script to estimate the dimension of the skein module of the
3-torus. Works by increasing the size of an (almost) triangular lattice shell:
all points in a box above the line y = -x, excluding a half-line to ensure this
is a subset of a fundamental domain for half-turn rotation. E.g for n = 4:

                        x x x x x x x x x
                          x x x x x x x x
                            x x x x x x x
                              x x x x x x
                                0 x x x x
                                    x x x
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
    Returns a tuple consisting of: a dictionary giving an order to the lattice
    points in a triangular shell, and a list giving the points in order.
    Ordering given by: (a, b) < (c, d) if a < c; (a, b) < (a, c) if b > c.
    That is, ordering goes left-right, top-bottom (lrtb) through the triangular
    shell.
    '''
    order_dict = {}
    points_in_order = []
    place = 0 # Place of the current lattice point in the order (to increment)

    # Incrementally loop over a, then decrement through b until we reach the
    # edge of the shell.
    for a in range(-1*shell_level, shell_level+1):
        b = shell_level
        if a <= 0:
            #Here the shell edge is the line y = -x; decrement until this point.
            while (a + b) >= 0:
                # Populate dict entry, increment place for next point,
                # decrement y coord.
                order_dict.update({(a, b): place})
                points_in_order.append(np.matrix([a, b]))
                place += 1
                b -= 1
        else:
            #Here the shell edge is the line y = -x + 1; decrement to this point
            while (a + b) >= 1:
                # Populate dict entry, increment place for next point,
                # decrement y coord.
                order_dict.update({(a, b): place})
                points_in_order.append(np.matrix([a, b]))
                place += 1
                b -= 1
    return (order_dict, points_in_order)

def get_relations(gamma, shell_level, order_func):
    '''
    Returns a list of linear relations between lattice points for a specified
    shell level. Each relation is a list representing a relation vector.

    Requires integer shell_level, and function order_func : int -> (dict, list)
    which should produce: a dictionary with keys being lattice points in a
    certain shell level, and values being their position in some sequential
    ordering specified by the function (this is required to store the relations
    as vectors); and a list of lattice points in order (required to pass back
    from a matrix entry to a lattice point easily, and to walk through the
    lattice).

    Performs a double loop through the lattice, obtains the relation for each
    pair of points, and discards trivial or out-of-range relations.
    '''

    relations = []
    N = (2*shell_level + 1)*(shell_level + 1) - shell_level # Total lattice pts.
    ordering = order_func(shell_level)[0] # Dictionary giving points an order.
    points_in_order = order_func(shell_level)[1] # Points listed in order.

    a = gamma[0, 0]
    b = gamma[0, 1]
    c = gamma[1, 0]
    d = gamma[1, 1]

    for p_0 in points_in_order:
        for p_1 in points_in_order:
            r = p_0[0, 0]
            s = p_0[0, 1]
            t = p_1[0, 0]
            u = p_1[0, 1]

            K = (-r*(r-1)*a*c - s*(s-1)*b*d)/2 - r*s*c*b # A constant

            # Rewrite indices for the generators appearing in the linear
            # relation coming from points ((a, b), (c, d)).
            x_0 = p_0 + p_1
            x_1 = p_0 - p_1
            x_2 = p_0 + np.matmul(p_1, gamma.T)
            x_3 = p_0 - np.matmul(p_1, gamma.T)

            x_0_tuple = tuple(x_0.tolist()[0])
            x_1_tuple = tuple(x_1.tolist()[0])
            x_2_tuple = tuple(x_2.tolist()[0])
            x_3_tuple = tuple(x_3.tolist()[0])

            # Check the relations are not out of range.
            if x_0_tuple in ordering.keys() and x_1_tuple in ordering.keys() and x_2_tuple in ordering.keys() and x_3_tuple in ordering.keys():
                #Create vectors corresponding to the four points
                x_0_vect = vector(QQ['q'].fraction_field(), [1 if i == ordering[x_0_tuple] else 0 for i in range(N)])
                x_1_vect = vector(QQ['q'].fraction_field(), [1 if i == ordering[x_1_tuple] else 0 for i in range(N)])
                x_2_vect = vector(QQ['q'].fraction_field(), [1 if i == ordering[x_2_tuple] else 0 for i in range(N)])
                x_3_vect = vector(QQ['q'].fraction_field(), [1 if i == ordering[x_3_tuple] else 0 for i in range(N)])

                # Compute the coefficients in the relations.
                Q_0 = q**(-s*t)
                Q_1 = q**(s*t)
                Q_2 = -q**(K - r*(c*t + d*u))
                Q_3 = -q**(K + r*(c*t + d*u))

                rel = Q_0*x_0_vect + Q_1*x_1_vect + Q_2*x_2_vect + Q_3*x_3_vect

                # Check the relations are not trivial.
                if not rel.is_zero():
                    relations.append(rel)
    return relations

import numpy as np
# Solicit user input
print("TWISTED TORUS SKEIN DIMENSION ESTIMATOR")
print("Input an SL_2(Z)-matrix to define a twisted 3-torus. Enter the matrix\n\n[[a b]\n [c d]]\n\nas the string a b c d and press return.")
user_input = [int(i) for i in input().split(" ")]
if len(user_input) != 4:
    print("Error! You did not enter 4 space-separated integers.")
    sys.exit(1)

gamma = np.matrix(user_input).reshape((2, 2))

if np.linalg.det(gamma) != 1:
    print("Error: the data you entered is not an SL_2(Z) matrix, must have determinant 1.")
    sys.exit(1)

print("You have entered the matrix \n\n[[%d %d]\n [%d %d]]\n" % tuple(user_input))

n = int(input("Enter the number of shell levels to check: "))#sage_eval(sys.argv[1]) # Number of levels of shell to compute for.
q = var('q') # Declare an indeterminate q.

for shell_level in range(n+1):
    # For each shell level, compute #{lattice points}.
    N = (2*shell_level + 1)*(shell_level + 1) - shell_level

    print("Calculating relations for level %d (%d lattice points) ..." % (shell_level, N))
    relations = get_relations(gamma, shell_level, order_lrtb)
    print("Found %d (non-independent) relations. Reducing ..." % len(relations))

    # Form a relation matrix, compute its rank; the dimension estimate is the
    # co-rank.
    A = matrix(QQ['q'].fraction_field(), relations)
    #print(A)
    relations_rank = A.rank()
    dim_estimate = N - relations_rank
    print("Dimension estimate for level %d: %d.\n" % (shell_level, dim_estimate))
