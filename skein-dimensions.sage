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
'''

import os
import sys
import numpy as np
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
                points_in_order.append(vector(QQ['q'].fraction_field(), [a, b]))
                place += 1
                b -= 1
        else:
            #Here the shell edge is the line y = -x + 1; decrement to this point
            while (a + b) >= 1:
                # Populate dict entry, increment place for next point,
                # decrement y coord.
                order_dict.update({(a, b): place})
                points_in_order.append(vector(QQ['q'].fraction_field(), [a, b]))
                place += 1
                b -= 1
    return (order_dict, points_in_order)

def get_relations_empty(gamma, shell_level, order_func):
    '''
    Returns a list of linear relations between lattice points for a specified
    shell level.

    Each ordered pair of lattice points determines a relation between four other
    lattice points, where lattice points correspond to generators of the skein
    module.

    Requires integer shell_level, and function order_func : int -> (dict, list)
    which should produce: a dictionary with keys being lattice points in a
    certain shell level, and values being their position in some sequential
    ordering specified by the function (this is required to map lattice points
    to indices of vectors in the space they span); and a list of lattice points
    in this order (required to produce the four related lattice points using
    basic linalg).

    Performs a double loop through the lattice, obtains the relation between
    four points for each pair of points, and discards trivial or out-of-range
    relations.
    '''

    relations = []
    N = (2*shell_level + 1)*(shell_level + 1) - shell_level # Total lattice pts.
    ordering, points_in_order = order_func(shell_level) # Dict and list of order

    #Unpack the matrix gamma.
    a = gamma[0, 0]
    b = gamma[0, 1]
    c = gamma[1, 0]
    d = gamma[1, 1]

    # Give a sage matrix version of gamma.
    Gamma = matrix(QQ['q'].fraction_field(), gamma.tolist())

    for p_0 in points_in_order:
        for p_1 in points_in_order:
            #Unpack the points
            r = p_0[0]
            s = p_0[1]
            t = p_1[0]
            u = p_1[1]

            #A constant appearing in our coefficients, we compute it in advance
            K = (-r*(r-1)*a*c - s*(s-1)*b*d)/2 - r*s*c*b

            # The linear relation is between the four lattice points below:
            x_0 = p_0 + p_1
            x_1 = p_0 - p_1
            x_2 = p_0 + p_1*Gamma.T
            x_3 = p_0 - p_1*Gamma.T

            #Get tuple versions of the above (to access the order dict)
            x_0_tuple = tuple([i for i in x_0])
            x_1_tuple = tuple([i for i in x_1])
            x_2_tuple = tuple([i for i in x_2])
            x_3_tuple = tuple([i for i in x_3])

            # Check the relations are not out of range.
            if x_0_tuple in ordering.keys() and x_1_tuple in ordering.keys() and x_2_tuple in ordering.keys() and x_3_tuple in ordering.keys():
                #Create vectors corresponding to the four lattice points.
                x_0_vect = vector(QQ['q'].fraction_field(), [1 if i == ordering[x_0_tuple] else 0 for i in range(N)])
                x_1_vect = vector(QQ['q'].fraction_field(), [1 if i == ordering[x_1_tuple] else 0 for i in range(N)])
                x_2_vect = vector(QQ['q'].fraction_field(), [1 if i == ordering[x_2_tuple] else 0 for i in range(N)])
                x_3_vect = vector(QQ['q'].fraction_field(), [1 if i == ordering[x_3_tuple] else 0 for i in range(N)])

                # Compute the coefficients in the relation.
                Q_0 = q**(-s*t)
                Q_1 = q**(s*t)
                Q_2 = -q**(K - r*(c*t + d*u))
                Q_3 = -q**(K + r*(c*t + d*u))

                #The relation is the following:
                rel = Q_0*x_0_vect + Q_1*x_1_vect + Q_2*x_2_vect + Q_3*x_3_vect

                # Check the relation is not trivial, then append.
                if not rel.is_zero():
                    relations.append(rel)
    return relations

def print_generators(shell_level, pivots, order_func):
    '''
    Prints a visualisation of the spanning lattice points to the command line.
    A spanning vector is denoted x in the lattice, other points are denoted .
    and axes are drawn using | and -. The visualisation is printed to the
    terminal (provided the terminal is large enough).
    Takes the shell level, a tuple giving the indices of the pivots of the
    relation matrix (these are the complement of the spanning set), and the
    order_func to map lattice points to indices for comparison.
    '''
    max_width = os.get_terminal_size().columns #Check terminal is wide enough
    if 2*(2*shell_level + 1) > max_width:
        print("Cannot display spanning set graphically.")
    else:
        ordering = order_func(shell_level)[0] # Dictionary giving points order.
        # Walk through (part of) the lattice Z^2 row by row, left to right.
        for y in range(shell_level, -1*shell_level - 1, -1):
            for x in range(-1*shell_level, shell_level + 1):
                # If a point is in the shell, check if it is NOT a pivot of the
                # relation matrix.
                if (x, y) in ordering.keys():
                    if not ordering[(x, y)] in pivots:
                        print("x ", end="") #Place an x for spanning vectors.
                    elif x == 0 and y == 0:
                        print("+ ", end="") #Origin.
                    elif x == 0:
                        print("| ", end="") #Y axis.
                    elif y == 0:
                        print("- ", end="") #X axis.
                    else:
                        print(". ", end="") #Generic lattice point.
                elif x == 0 and y == 0:
                    print("+ ", end="") #Origin.
                elif x == 0:
                    print("| ", end="") #Y axis.
                elif y == 0:
                    print("- ", end="") #X axis.
                else:
                    print(". ", end="") #Generic lattice point.
            print("") # Complete line with \n.
        print("") # Pad below.
    return None

def order_lexi():
    Z_2 = Integers(2)
    order_dict = {}
    in_order = []
    place = 0
    for a_0 in Z_2:
        for a_1 in Z_2:
            order_dict.update({(a_0, a_1) : place})
            in_order.append(vector([a_0, a_1]))
            place += 1

    return (order_dict, in_order)

def get_dim_single_skein(gamma):
    '''
    Takes a 4x4 np array gamma, and returns the dimension of the single skein
    part of the twisted torus, where gamma is an SL_2(Z)-matrix defining the
    twisting.

    The vector space in question is a quotient of C[X, Y]/(X^2, Y^2) by
    some relations. The basis {1, X, Y, XY} is represented by vectors of length
    2 with Z/2-entries, i.e. X^aY^b is [a, b], and these are ordered
    lexicographically.

    The implementation is similar to get_relations.
    '''
    Z_2 = Integers(2) # Integers mod 2
    order_dict, in_order = order_lexi() # A dict and list to order the basis
    Gamma = matrix(Z_2, gamma.tolist()) # mod 2 sage version of gamma

    rels = [] # Ready to record the twisted commutator relations.

    # Walk thru pairs of basis elements (R, S) viewed as 1x2 Z/2 vectors (r, s)
    for r in in_order:
        for s in in_order:
            first_term = r + s # First term in the commutator is RS, i.e. r + s
            second_term = s*Gamma.T + r # Second term is twisted by Gamma

            # Turn these into tuples to access the order dict
            first_term_tuple = tuple([i for i in first_term])
            second_term_tuple = tuple([i for i in second_term])

            # Use the order dict to view each term of the relation as a vector
            # in a 4d vecctor space.
            first_term_vector = vector(QQ, [1 if i == order_dict[first_term_tuple] else 0 for i in range(4)])
            second_term_vector = vector(QQ, [1 if i == order_dict[second_term_tuple] else 0 for i in range(4)])

            #Get the relation, and append it if non-trivial.
            rel = first_term_vector - second_term_vector
            if not rel.is_zero():
                rels.append(rel)

    # Calclulate the rank of the relation matrix: the dimention of the single
    # skein part is the co-rank.
    R = matrix(QQ, rels)
    dim = 4 - R.rank()

    return dim

# Solicit user input
print("TWISTED TORUS SKEIN DIMENSION ESTIMATOR")
print("Input an SL_2(Z)-matrix to define a twisted 3-torus. Enter the matrix\n\n[[a b]\n [c d]]\n\nas the string a b c d and press return.")
user_input = [int(i) for i in input().split(" ")]
if len(user_input) != 4:
    print("Error! You did not enter 4 space-separated integers.")
    sys.exit(1)

# This matrix defines the twisted torus
gamma = np.matrix(user_input).reshape((2, 2))
print("You have entered the matrix \n\n[[%d %d]\n [%d %d]]\n" % tuple(user_input))

if np.linalg.det(gamma) != 1:
    print("Error: the data you entered is not an SL_2(Z) matrix, must have determinant 1.")
    sys.exit(1)

n = int(input("Enter the number of shell levels to check: "))

# Declare an indeterminate q.
q = var('q')

# Estimate the skein module dimension for each shell level.
for shell_level in range(n+1):
    # For each shell level, compute #{lattice points}.
    N = (2*shell_level + 1)*(shell_level + 1) - shell_level
    print("Calculating relations for level %d (%d lattice points) ..." % (shell_level, N))
    relations = get_relations_empty(gamma, shell_level, order_lrtb)
    print("Found %d (non-independent) relations. Reducing ..." % len(relations))
    # Form a relation matrix, compute its pivots; the dimension estimate is the
    # co-rank.
    A = matrix(QQ['q'].fraction_field(), relations)
    pivots = A.pivots()
    dim_estimate = N - len(pivots)
    print("Dimension estimate for level %d: %d.\n\nVisualisation:\n" % (shell_level, dim_estimate))
    print_generators(shell_level, pivots, order_lrtb)
