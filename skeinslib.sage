'''
USAGE In a sage interactive session or a sage script, load these functions using

load("skeinslib.sage")

*Requirements* Sage v9.2 or higher; Python v3.7 or higher.

*Note* that the function write_dim_table required pandas imported as pd. This
can be installed by opening a Sage shell and running

pip install pandas

or by running

sage -pip install pandas

at the command line.

DESCRIPTION A library of functions used to estimate the dimension of the skein
module of the twisted 3-torus. Gives the dimension of the single skein part, and
estimates the dimension of the empty skein part.

The empty skein estimate works by increasing the size of an (almost) triangular
lattice shell: all points in a box above the line y = -x, excluding a half-line
to ensure this is a subset of a fundamental domain for half-turn rotation.
E.g for n = 4:

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
gives an upper bound on the dimension of the empty skein part of the skein
module.

The dimension of the single skein part is easily computed, and there is no
estimation in this figure.
'''

import csv
import itertools
import os
import pandas as pd

from sage.all import *

var('q')
K = FractionField(PolynomialRing(QQ, 'q', sparse=True))

def order_by_shell_level(shell_level):
    '''
    Returns a dictionary giving an order to the lattice points in a triangular
    shell. The ordering within the dictionary (dict order persists from Python
    3.7) is the ordering, and the dict is used for a fast lookup of ordering of
    the points.
    Ordering is given by shell level, and within this goes left-right,
    top-bottom. E.g

                    13 14 15 16 17 18 19
                        5  6  7  8  9 20
                           1  2  3 10 21
                              0  4 11 22
                                   12 23
                                      24


    The order is to use the points to index a basis of the space which is
    quotiented by the twisted commutator relations. The keys in the order dict
    are immutable sage vectors in ZZ, values are position in ordering.
    Implemented recursively.
    '''

    # Base case: shell level 0
    if shell_level == 0:
        order_dict = {vector(ZZ, [0, 0], immutable=True) : 0}

    # Otherwise, recurse
    else:
        order_dict = order_by_shell_level(shell_level - 1)
        place = len(order_dict.keys())

        # We will walk through points (a, b) in the L shape of points in shell
        # level n but not level n-1.
        b = shell_level
        a = -1*shell_level

        # We turn the corner for a = b = shell_level.
        while a < shell_level:
             order_dict.update({vector(ZZ, [a, b], immutable=True) : place})
             place += 1
             a += 1

        #Here the shell edge is the line y = -x + 1; decrement to this point
        while b > -1*shell_level:
            order_dict.update({vector(ZZ, [a, b], immutable=True) : place})
            place += 1
            b -= 1

    return order_dict

def get_relations_empty(gamma, shell_level, order_func):
    '''
    Returns a dict of linear relations between lattice points for a specified
    shell level.

    Each ordered pair of lattice points determines a relation between four other
    lattice points, where lattice points correspond to generators of the empty
    part of the skein module.

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

    relations = {}
    N = (2*shell_level + 1)*(shell_level + 1) - shell_level # Total lattice pts.
    ordering = order_func(shell_level) # Dict and list of order

    #Unpack the matrix gamma.
    a = gamma[0, 0]
    b = gamma[0, 1]
    c = gamma[1, 0]
    d = gamma[1, 1]

    q = var('q') # Must be defined here to alllow compiled sage.

    for p_0 in ordering.keys():
        for p_1 in ordering.keys():
            #Unpack the points
            r = p_0[0]
            s = p_0[1]
            t = p_1[0]
            u = p_1[1]

            #A constant appearing in our coefficients, we compute it in advance
            C = (-t*(t-1)*a*c - u*(u-1)*b*d)/2 - t*u*c*b

            # The linear relation is between the four lattice points below:
            x_0 = vector(ZZ, p_0 + p_1, immutable=True)
            x_1 = vector(ZZ, p_0 - p_1, immutable=True)
            x_2 = vector(ZZ, p_0 + p_1*gamma.T, immutable=True)
            x_3 = vector(ZZ, p_0 - p_1*gamma.T, immutable=True)

            # Check the relations are not out of range.
            if x_0 in ordering.keys() and x_1 in ordering.keys() and x_2 in ordering.keys() and x_3 in ordering.keys():
                #Create vectors corresponding to the four lattice points.
                x_0_vect = vector(K, [1 if i == ordering[x_0] else 0 for i in range(N)], sparse=True)
                x_1_vect = vector(K, [1 if i == ordering[x_1] else 0 for i in range(N)], sparse=True)
                x_2_vect = vector(K, [1 if i == ordering[x_2] else 0 for i in range(N)], sparse=True)
                x_3_vect = vector(K, [1 if i == ordering[x_3] else 0 for i in range(N)], sparse=True)

                # Compute the coefficients in the relation.
                Q_0 = q**(-s*t)
                Q_1 = q**(s*t)
                Q_2 = -q**(C - r*(c*t + d*u))
                Q_3 = -q**(C + r*(c*t + d*u))

                #The relation is the following:
                rel = Q_0*x_0_vect + Q_1*x_1_vect + Q_2*x_2_vect + Q_3*x_3_vect

                # Check the relation is not trivial, then append.
                if not rel.is_zero():
                    relations[(p_0, p_1)] = rel

    return relations

def get_new_relations_empty(gamma, shell_level, order_func):
    '''
    Returns a list of linear relations between lattice points for a specified
    shell level. Like get_relations_empty except only returns the relations
    which are new at the specified shell level (did not occur at the previous
    level).

    Each ordered pair of lattice points determines a relation between four other
    lattice points, where lattice points correspond to generators of the empty
    part of the skein module.

    Requires integer shell_level, and function order_func : int -> (dict, list)
    which should produce: a dictionary with keys being lattice points in a
    certain shell level, and values being their position in some sequential
    ordering specified by the function (this is required to map lattice points
    to indices of vectors in the space they span); and a list of lattice points
    in this order (required to produce the four related lattice points using
    basic linalg).

    Performs a double loop through the lattice, obtains the relation between
    four points for each pair of points, and discards trivial or out-of-range
    relations, or relations which lie in the previous shell level.
    '''

    relations = []
    N = (2*shell_level + 1)*(shell_level + 1) - shell_level # Total lattice pts.
    ordering = order_func(shell_level) # Dict and list of order
    M = (2*shell_level - 1)*(shell_level) - shell_level + 1 # First elements from the previous shell level

    #Unpack the matrix gamma.
    a = gamma[0, 0]
    b = gamma[0, 1]
    c = gamma[1, 0]
    d = gamma[1, 1]

    q = var('q') # Must be defined here to alllow compiled sage.

    for p_0 in ordering.keys():
        for p_1 in ordering.keys():
            #Unpack the points
            r = p_0[0]
            s = p_0[1]
            t = p_1[0]
            u = p_1[1]

            #A constant appearing in our coefficients, we compute it in advance
            C = (-t*(t-1)*a*c - u*(u-1)*b*d)/2 - t*u*c*b

            # The linear relation is between the four lattice points below:
            x_0 = vector(ZZ, p_0 + p_1, immutable=True)
            x_1 = vector(ZZ, p_0 - p_1, immutable=True)
            x_2 = vector(ZZ, p_0 + p_1*gamma.T, immutable=True)
            x_3 = vector(ZZ, p_0 - p_1*gamma.T, immutable=True)

            # Check the relations are not out of range.

            if x_0 in ordering and x_1 in ordering and x_2 in ordering and x_3 in ordering:
                if not (x_0  in list(ordering.keys())[:M] and x_1 in list(ordering.keys())[:M] and x_2 in list(ordering.keys())[:M] and x_3 in list(ordering.keys())[:M]):

                    x_0_vect = vector(K, [1 if i == ordering[x_0] else 0 for i in range(N)], sparse=True)
                    x_1_vect = vector(K, [1 if i == ordering[x_1] else 0 for i in range(N)], sparse=True)
                    x_2_vect = vector(K, [1 if i == ordering[x_2] else 0 for i in range(N)], sparse=True)
                    x_3_vect = vector(K, [1 if i == ordering[x_3] else 0 for i in range(N)], sparse=True)

                    # Compute the coefficients in the relation.
                    Q_0 = q**(-s*t)
                    Q_1 = q**(s*t)
                    Q_2 = -q**(C - r*(c*t + d*u))
                    Q_3 = -q**(C + r*(c*t + d*u))

                    #The relation is the following:
                    rel = Q_0*x_0_vect + Q_1*x_1_vect + Q_2*x_2_vect + Q_3*x_3_vect

                    # Check the relation is not trivial, then append.
                    if not rel.is_zero():
                        relations.append(rel)

    return relations

def print_generators(shell_level, spanning_set, order_func):
    '''
    Prints a visualisation of the spanning lattice points to the command line.
    A spanning vector is denoted x in the lattice, other points are denoted .
    and axes are drawn using | and -. The visualisation is printed to the
    terminal (provided the terminal is large enough).
    Takes the shell level, a tuple giving the indices of the spanning set, and
    the order_func to map lattice points to indices for comparison.
    '''

    max_width = os.get_terminal_size().columns #Check terminal is wide enough
    if 2*(2*shell_level + 1) > max_width:
        print("Cannot display spanning set graphically.")
    else:
        ordering = order_func(shell_level) # Dictionary giving points order.
        # Walk through (part of) the lattice Z^2 row by row, left to right.
        for y in range(shell_level, -1*shell_level - 1, -1):
            for x in range(-1*shell_level, shell_level + 1):
                # If a point is in the shell, check if it is NOT a pivot of the
                # relation matrix.
                if vector(ZZ, [x, y], immutable=True) in ordering.keys():
                    if ordering[vector(ZZ, [x, y], immutable=True)] in spanning_set:
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
    '''
    Returns an ordering on 1x2 vector representations of basis elements of
    C[X, Y]/(X^2 - 1, Y^2 - 1). The element R = X^aY^b is the vector r = [a, b],
    so that RS is given by rs and gamma.R is r*gamma.T, and the basis is 
    {1, X, Y, XY}.

    Returns a dictionary giving lexicographical ordering on these vectors.

    order_dict: keys are 2d vectors (with Z/2 entries), values are their place
    in the lexicographical ordering: used to pass from r to the 4d vector R. The
    keys must be immutable to be hashable. The order dict is itself in order,
    and the ordering persists (from Python 3.7), i.e. oder-dict.keys() is an
    iterable with the vectors given in the ordering.
    '''

    Z_2 = Integers(2) # Work mod 2
    order_dict = {}
    place = 0 # Record the place in the ordering

    # Iterate over (Z/2)^2 in lexicographical order
    for a_0 in Z_2:
        for a_1 in Z_2:
            # Update the dict and list
            order_dict.update({vector(Z_2, [a_0, a_1], immutable=True) : place})
            place += 1

    return order_dict

def get_dim_single_skein(gamma):
    '''
    Takes an SL_2(Z) matrix gamma, and returns the dimension of the single skein
    part of the twisted torus, where gamma is an SL_2(Z)-matrix defining the
    twisting.

    The vector space in question is a quotient of C[X, Y]/(X^2 - 1 , Y^2 - 1) by
    some relations. The basis {1, X, Y, XY} is represented by vectors of length
    2 with Z/2-entries, i.e. X^aY^b is [a, b], and these are ordered
    lexicographically.

    The implementation is similar to get_relations_empty: for all pairs of basis
    elements of C[X, Y]/(X^2 - 1, Y^2 - 1) we get the gamma-twisted commutators, 
    then the corank of these relations is the required dimension.
    '''

    Z_2 = Integers(2) # Integers mod 2
    order_dict = order_lexi() # A dict and list to order the basis

    rels = [] # Ready to record the twisted commutator relations.

    # Walk thru pairs of basis elements (R, S) viewed as 1x2 Z/2 vectors (r, s)
    for r in order_dict.keys():
        for s in order_dict.keys():
            first_term = vector(Z_2, r + s, immutable=True) # First term in the commutator is RS, i.e. r + s
            second_term = vector(Z_2, s*gamma.T + r, immutable=True) # Second term is twisted by Gamma

            # Use the order dict to view each term of the relation as a vector
            # in a 4d vecctor space.
            first_term_vector = vector(QQ, [1 if i == order_dict[first_term] else 0 for i in range(4)])
            second_term_vector = vector(QQ, [1 if i == order_dict[second_term] else 0 for i in range(4)])

            #Get the relation, and append it if non-trivial.
            rel = first_term_vector - second_term_vector
            if not rel.is_zero():
                rels.append(rel)

    # Calclulate the rank of the relation matrix: the dimention of the single
    # skein part is the co-rank.
    R = matrix(QQ, rels)
    dim = 4 - R.rank()

    return dim

def compute_reduced_matrix_defrost(gamma, shell_level, data_path):
    '''
    Invokes compute_reduced_matrix, but only computing to a specified recursion
    base which is stored in the .sobj file given by data_path. This must unpack
    to a list/tuple of length 3, with contents [shell_levels, A, A_reduced, [dim_estimates]] from
    a previous run of compute_reduced_matrix.

    Note that, if the default file naming conventions of the program are used,
    then the file at data_path will be overwritten when this function is called.
    '''

    if not os.exists(data_path):
        print("Error, file not found")
        sys.exit(1)

    data = sage.misc.persist.load("data_path")

    return compute_reduced_matrix(gamma, shell_level, interactive_flag=False, base_level = data[0], base_data = data[1:])


def compute_reduced_matrix(gamma, shell_level, interactive_flag, base_level = 1, base_data=None):
    '''
    Computes the reduced matrix of relations at a given shell level, as well as
    its co-rank, which is an estimate of the dimension of the empty part of the
    skein module

    Returns a tuple (A, A_reduced, dim_estimates), where A is the matrix and
    dim_estimates is the list of estimated dimensions for shell levels <= to the
    shell level parameter.

    Works recursively: the relation matrix for level n is built up as

                [[  X  | 0 ]
                 [---------]
                 [relations]]

    where X is the matrix of relations for level n-1, in REF, obtained from a
    recursive call to the function, and the list relations is a list of
    relations at level n that were not already included at level n-1.

    Supports recursion to a specified base_level, with base_data the results of
    a previous call to the function (that is, a tuple of length 3). Usually this
    would be loaded as a sage object saved from a previous session. If the
    recursion base is unspecified then the default is recurse to shell_level 1
    and compute this data from scratch.
    '''

    # Base case.
    if shell_level == base_level:
        if base_level == 1: # Recursing all the way to level 1
            # For each shell level, compute #{lattice points}.
            N = (2*shell_level + 1)*(shell_level + 1) - shell_level
            # And the number of points in previous level.
            M = (2*shell_level - 1)*(shell_level) - shell_level + 1
            if interactive_flag:
                print("Calculating relations for level %d (%d lattice points) ..." % (shell_level, N))

            # Get the relations new to this shell level.
            relations = get_new_relations_empty(gamma, shell_level, order_by_shell_level)
            if interactive_flag:
                print("Found %d (non-independent) relations.\n" % len(relations))
            # Form a relation matrix, and reduce it.
            A = matrix(K, relations, sparse=True, immutable=True)
            A_reduced = A.rref()

            # Get the spanning set
            ordering = order_by_shell_level(shell_level)
            spanning_set = get_spanning_set(A_reduced, ordering, shell_level)

            # Base of the recursion: return the reduced matrix and spanning set.
            dim_estimate = len(spanning_set)
            if interactive_flag:
                print("Dimension estimate for empty skein part at level %d: %d.\n\nVisualisation:\n" % (shell_level, dim_estimate))
                print_generators(shell_level, spanning_set, order_by_shell_level)

            return (A, A_reduced, [dim_estimate])

        else: # Recursing to a higher level and picking up where we left off.
            if interactive_flag:
                print("Using previosuly calculated relations for level %d..." % (shell_level))
            return base_data

    else:
        # For each shell level, compute #{lattice points}.
        N = (2*shell_level + 1)*(shell_level + 1) - shell_level
        # And the number of points in previous level.
        M = (2*shell_level - 1)*(shell_level) - shell_level + 1
        if interactive_flag:
            print("Calculating relations for level %d (%d lattice points) ..." % (shell_level, N))

        # Get the relations new to this shell level.
        relations = get_new_relations_empty(gamma, shell_level, order_by_shell_level)
        if interactive_flag:
            print("Found %d (non-independent) relations. Reducing ..." % len(relations))
        # Get the relation matrix for the previous shell level, recursively.
        prev_data = compute_reduced_matrix(gamma, shell_level - 1, interactive_flag, base_level, base_data)
        A_old = prev_data[1]
        dimensions = prev_data[2]

        # Handle the case when there were no relations in the previous level
        if len(A_old.rows()) == 0:
            A = matrix(K, relations, sparse=True, immutable=True)
        else:
            # Use the new relations and old, reduced matrix to build the relations
            # matrix for this shell level, and reduce.
            Zeros_right = zero_matrix(K, A_old.nrows(), N - M, sparse=True)
            A_upper = block_matrix([[A_old, Zeros_right]])
            A_lower = matrix(K, relations, sparse=True, immutable=True)
            A = block_matrix([[A_upper], [A_lower]])

        A_reduced = A.rref()

        # Get the spanning set
        ordering = order_by_shell_level(shell_level)
        spanning_set = get_spanning_set(A_reduced, ordering, shell_level)

        # Add the dimension estimate to the list, print output, and return
        dim_estimate = len(spanning_set)
        dimensions.append(dim_estimate)
        if interactive_flag:
            print("Dimension estimate for empty skein part at level %d: %d.\n\nVisualisation:\n" % (shell_level, dim_estimate))
            print_generators(shell_level, spanning_set, order_by_shell_level)

        return (A, A_reduced, dimensions)


def get_spanning_set(A, ordering, shell_level):
    '''
    Given a relation matrix A, the ordering, and the shell level, returns a list
    of the indices of spanning vectors for the approximation of the empty skein
    part at this shell level.

    For each lattice vector in the ordering, checks if it is in the span of the
    matrix B given by the relations A and the previous vectors in the ordering.
    This is checked by augmenting B with the corresponding length N basis
    vector and comparing with the rank of B. If the lattice point corresponds to
    a vector not in the span of B, its position in the ordering is added to the
    list of spanning vectors.
    '''

    N = (2*shell_level + 1)*(shell_level + 1) - shell_level
    spanning_set = []

    #Handle the case of the empty relation matrix
    if len(A.rows()) == 0:
        spanning_set = list(ordering.values())

    else:
        # Iterate through the lattice.
        for lattice_pt in ordering.keys():
            #Corresponding length N vector.
            basis_elt = vector(K, [1 if i == ordering[lattice_pt] else 0 for i in range(N)], sparse=True)
            A = A.rref() # Reduce A at each step.
            rk_A = A.rank() # Store the rank of A
            A = block_matrix([[A], [matrix(basis_elt)]]) # Augment A with the new vector.
            # If not in this span, add the position of the lattice point to the list
            # of indices of spanning vectors (used to print generators).
            if A.rank() > rk_A:
                spanning_set.append(ordering[lattice_pt])

    return spanning_set

def compute_and_write(sequence, gamma, shell_levels, try_defrost, dir_in, dir_out, output_path, cache_path):
    '''
    A helper function for generate_raw_data, handles the subroutine of
    collating the results of dimension computations for M (up to shell_levels
    cutoff) and writes to the file at output_path, in append mode.

    Also stores the relation matrix, reduced version, and dim estimates as a 
    tuple sage object, and maintains a persistent cache. If try_defrost=True, 
    will attempt to load the results of previous computations and use these 
    to only recurse to a specified base in compute-reduced_matrix.
    '''

    # Store the computed sage objects in dir_out, indexed by the matrix entries.
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)

    seq_string = "{0}_{1}_{2}_{3}".format(gamma[0, 0], gamma[0, 1], gamma[1, 0], gamma[1, 1])
    seq_path_in = os.path.join(dir_in, seq_string+".sobj")
    seq_path_out = os.path.join(dir_out, seq_string+".sobj")

    # Declare an indeterminate q.
    q = var('q')

    # Get the dimension of the single skein.
    dim_single = get_dim_single_skein(gamma)

    # If we are trying to defrost computations, and the data exists, use this to continue our computations
    # of the data for the empty skein part. Otherwise recurse to level 0.
    if try_defrost and os.path.exists(seq_path_in):
        data = sage.misc.persist.load(seq_path_in)
        empty_data = compute_reduced_matrix(gamma, shell_levels, False, base_level=data[0], base_data=data[1:])
    else:
        empty_data = compute_reduced_matrix(gamma, shell_levels, False)

    # Extract some of the numbers we are interested in.
    dim_estimates = empty_data[2]
    dim_total = dim_single + dim_estimates[-1]

    # Write the relevant data to the output files.
    with open(output_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([gamma.trace(), gamma[0, 0], gamma[0, 1], gamma[1, 0], gamma[1, 1], dim_single, dim_total] + dim_estimates + list(sequence))
    f.close()

    A = empty_data[0]
    A_reduced = empty_data[1]
    sage.misc.persist.save([shell_levels, A, A_reduced, dim_estimates], seq_path_out)

    # Write the computed sequence to the persistent cache file.
    with open(cache_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(list(sequence))

    return None

def compute_write_low_trace(shell_levels, try_defrost, dir_in, dir_out, output_path, cache_path):
    '''
    Perform the computations for matrices of trace with abs val < 2, with a
    specified max shell level, and write to the file at path.
    '''

    #Trace 0 matrix
    S = matrix(ZZ, 2, [0, -1, 1, 0])

    #Trace 1 matrices
    M_0 = matrix(ZZ, 2, [1, -1, 1, 0])
    M_1 = matrix(ZZ, 2, [0, 1, -1, 1])

    low_trace = [S, M_0, M_1]

    # Compute dimensions for the 3 matrices of low trace.
    for M in low_trace:
        compute_and_write([], M, shell_levels, try_defrost=try_defrost, dir_in=dir_in, dir_out=dir_out, output_path=output_path, cache_path=cache_path)

    return None

def compute_write_from_seq(sequence, shell_levels, try_defrost, dir_in, dir_out, output_path, cache_path):
    '''
    Performs the dimension computations for an SL_2(Z) matrix of form
        R^{a_1}L^{a_2}...(R or L)^{a_k}
    given a sequence (a_k). Computations performed up to a max shell level, and
    written to path.

    Here R = [[1, 1], [0, 1]], L = [[1, 0], [1, 1]] so that R^n is
    [[1, n], [0, 1]] and this is how we implement the exponentiation (similar
    for L).
    '''

    #Generators for |trace| >= 2 matrices
    R = matrix(ZZ, 2, [1, 1, 0, 1])
    L = matrix(ZZ, 2, [1, 0, 1, 1])

    M = matrix(ZZ, 2, [1, 0, 0, 1])
    # Build up the word indexed by this sequence
    for i in range(len(sequence)):
        if i % 2 == 0:
            #Multiply by R^a
            M = M*matrix(ZZ, 2, [1, sequence[i], 0, 1])
        else:
            #Multiply by L^a
            M = M*matrix(ZZ, 2, [1, 0, sequence[i], 1])

    compute_and_write(sequence=sequence, gamma=M, shell_levels=shell_levels, try_defrost=try_defrost, dir_in=dir_in, dir_out=dir_out, output_path=output_path, cache_path=cache_path)
    return None

def seq_has_been_checked(seq, cache):
    '''
    Checks is the given sequence, or any rotation thereof, is in the cache of
    already-checked sequences.
    '''

    for i in range(len(seq)):
        seq = seq[1:] + seq[:1] # Rotate the sequence once.
        if seq in cache: # Check cache membership.
            return True
    return False

def generate_raw_data(shell_levels, append=False, try_defrost=True, dir_in="data", dir_out="data", output_path="skeindims-rawdata.csv", cache_path="seq-cache.csv"):
    '''
    Generate several SL_2(Z) matrices, compute their skein dimension estimates,
    and write this data to a csv file.

    Maintains a cache of previously checked sequences. In append mode, this is
    loaded from a csv and the new dimension data is appended to an existing
    output file.

    All data is stored in the directory specified by dir_out. 
    
    In try_defrost mode, this dir_in will be searched for .sobj files containing
    the results of previous computations, in which case the computation will recurse
    to this base step.

    Shell_levels is the total number of shell levels, including level 0 (!)
    '''

    half_max_seq_len = 2

    if not append:
        # Open a file for the raw data, and write a header.
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["trace", "a", "b", "c", "d", "single_dim", "total_dim"] + ["shell_{n}".format(n=i) for i in range(shell_levels)] + ["seq_{n}".format(n=i) for i in range(2*half_max_seq_len)])
        # Reset the cache
        if os.path.exists(cache_path):
            os.remove(cache_path)

    # Dimensions for low trace matrices.
    compute_write_low_trace(shell_levels=shell_levels, try_defrost=try_defrost, dir_in=dir_in, dir_out=dir_out, output_path=output_path, cache_path=cache_path)

    # Dimensions for the family of shears (|trace = 2|).
    for n in range(11):
        compute_write_from_seq(sequence=[n], shell_levels=shell_levels, try_defrost=try_defrost, dir_in=dir_in, dir_out=dir_out, output_path=output_path, cache_path=cache_path)

    # Dimensions for matrices of |trace| >  2

    
    cache = []  # Previously checked sequences.
    if append:
        if os.path.exists(cache_path):
            with open(cache_path, newline="") as cache_file:
                cache_reader = csv.reader(cache_file)
                cache = [tuple([int(s) for s in seq]) for seq in cache_reader]

    # The space of sequences to search.
    for half_seq_len in range(1, half_max_seq_len + 1):
        seq_len = 2*half_seq_len # Sequence lengths must be even.
        sequences = itertools.product(range(1, 11), repeat=seq_len)
        for sequence in sequences:
            #Exclude previously checked sequences up to cyclic permutation.
            if not seq_has_been_checked(sequence, cache):
                compute_write_from_seq(sequence=sequence, shell_levels=shell_levels, try_defrost=try_defrost, dir_in=dir_in, dir_out=dir_out, output_path=output_path, cache_path=cache_path)
                cache.append(sequence)
    

    return None

def write_dim_table(rawpath, outpath, shell_levels):
    '''
    Parses data produced by generate_raw_data, stored in the file rawpath, and
    gives this as a formatted table written to the file at outpath.
    '''

    df  = pd.read_csv(rawpath, dtype=str)

    #Calculate values that may be required for padding and alignment.
    max_tr_len = df["trace"].map(len).max()
    max_entry_len = df[["a", "b", "c", "d"]].applymap(len).values.max()
    max_empty_dim_len = df["shell_{n}".format(n=shell_levels-1)].map(len).max()
    max_total_dim_len = df["total_dim"].map(len).max()

    mat_line_length = 5 + max_entry_len*2

    with open(outpath, "w") as f:
        #Table header
        print("TRACE\t\tMATRIX" + " "*(mat_line_length-len("MATRIX")) + "\t\tSINGLE\tEMPTY\tTOTAL\t\t\tSHELL ESTIMATES", file=f)

        #Iterate through the entries and print rows.
        for idx, row in df.iterrows():
            print("{tr: >{max_tr_len}s}\t\t".format(tr=row["trace"], max_tr_len=max(max_tr_len, len("TRACE"))), end="", file=f)
            print("[[{a: >{max_entry_len}s} {b: >{max_entry_len}s}] \t\t".format(a=row["a"], b=row["b"], max_entry_len=max_entry_len), end="", file=f)
            print("{sing: >6s}\t".format(sing=row["single_dim"]), end="", file=f)
            print("{empty: >{max_empty_dim_len}s}\t".format(empty=row["shell_{n}".format(n=shell_levels-1)], max_empty_dim_len = max(max_empty_dim_len, len("EMPTY"))), end="", file=f)
            print("{tot: >{max_total_dim_len}s}\t\t\t".format(tot=row["total_dim"], max_total_dim_len=max(max_total_dim_len, len("TOTAL"))), end="", file=f)
            print("{ests}".format(ests=row[["shell_{n}".format(n=i) for i in range(shell_levels)]].values.tolist()), file=f)
            print(" "*max(max_tr_len,  5) + "\t\t" + " ", end="", file=f)
            print("[{c: >{max_entry_len}s} {d: >{max_entry_len}s}]]".format(c=row["c"], d=row["d"], max_entry_len=max_entry_len), file=f)
            print("", file=f)
        f.close()
