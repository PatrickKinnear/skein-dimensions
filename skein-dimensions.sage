#!/usr/bin/env sage

'''
USAGE Ensure SAGE_ROOT is stored in your PATH, and run

./skein-dimensions.sage mode [rawpath [outpath [shell-level]]]

where mode is a string specifying the mode, one of:

(i) interactive
(g) generation
(w) write
(gw) generate-and-write

*Note* that w and gw modes require pandas. This can be installed by opening a
Sage shell and running

pip install pandas

*Interactive Mode*
You will be prompted for an SL_2(Z)-matrix and an integer n, where n gives the
number of levels to search through in a triangular shell (see description), to
estimate the dimension of the empty skein part of the skein module of T^2 x S^1
twisted by the specified matrix. You may also give the letters I, S, or T to
compute for these matrices, or the same preceded by - (for the negative version)

*Generation Mode*
The program will generate matrices in SL_2(Z), grouped by trace, and compute the
dimension estimate of the skein module of the twisted torus defined by this
matrix. The data is written in raw form to a csv file in the working directory.
The output filepath may be specified at the command line. Default is
skein-dims-rawdata.csv. Can also specify shell level (defaults to 6).

*Write Mode*
The program will produce a table giving some SL_2(Z) matrices, the dimension of
the single skein part of the skein module of the twisted torus defined by this
matrix, and an estimate of the dimension of the empty skein part. The raw data
for the table is the output of a run of generation mode, so the program must be
run once in g mode before w mode is used. The path for the raw data file and the
formatted output file may be specified at the command line. Defaults are
skein-dims-rawdata.csv and skein-dims-printed.txt respectively. Can also specify
shell level (defaults to 6).

*Generate-Write Mode*
The result of running generation mode followed immediately by write mode.

DESCRIPTION A script to estimate the dimension of the skein module of the
twisted 3-torus. Gives the dimension of the single skein part, and estimates the
dimension of the empty skein part.

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

import os
import sys
import csv
import random as rand
import pandas as pd
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
                order_dict.update({vector(ZZ, [a, b], immutable=True): place})
                place += 1
                b -= 1
        else:
            #Here the shell edge is the line y = -x + 1; decrement to this point
            while (a + b) >= 1:
                # Populate dict entry, increment place for next point,
                # decrement y coord.
                order_dict.update({vector(ZZ, [a, b], immutable=True): place})
                place += 1
                b -= 1
    return order_dict

def get_relations_empty(gamma, shell_level, order_func):
    '''
    Returns a list of linear relations between lattice points for a specified
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

    relations = []
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
            K = (-r*(r-1)*a*c - s*(s-1)*b*d)/2 - r*s*c*b

            # The linear relation is between the four lattice points below:
            x_0 = vector(ZZ, p_0 + p_1, immutable=True)
            x_1 = vector(ZZ, p_0 - p_1, immutable=True)
            x_2 = vector(ZZ, p_0 + p_1*gamma.T, immutable=True)
            x_3 = vector(ZZ, p_0 - p_1*gamma.T, immutable=True)

            # Check the relations are not out of range.
            if x_0 in ordering.keys() and x_1 in ordering.keys() and x_2 in ordering.keys() and x_3 in ordering.keys():
                #Create vectors corresponding to the four lattice points.
                x_0_vect = vector(FractionField(PolynomialRing(QQ, 'q', sparse=True)), [1 if i == ordering[x_0] else 0 for i in range(N)], sparse=True)
                x_1_vect = vector(FractionField(PolynomialRing(QQ, 'q', sparse=True)), [1 if i == ordering[x_1] else 0 for i in range(N)], sparse=True)
                x_2_vect = vector(FractionField(PolynomialRing(QQ, 'q', sparse=True)), [1 if i == ordering[x_2] else 0 for i in range(N)], sparse=True)
                x_3_vect = vector(FractionField(PolynomialRing(QQ, 'q', sparse=True)), [1 if i == ordering[x_3] else 0 for i in range(N)], sparse=True)

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
        ordering = order_func(shell_level) # Dictionary giving points order.
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
    '''
    Returns an ordering on 1x2 vector representations of basis elements of
    C[X, Y]/(X^2, Y^2). The element R = X^aY^b is the vector r = [a, b], so that
    RS is given by rs and gamma.R is r*gamma.T, and the basis is {1, X, Y, XY}.

    Returns a tuple (order_dict, in_order) giving lexicographical ordering on
    these vectors.

    order_dict: keys are 2d vectors (with Z/2 entries), values are their place
    in the lexicographical ordering: used to pass from r to the 4d vector R. The
    keys must be tuples to be hashable.

    in_order: a list of the vectors, as 2d sage vectors, in lexicographical
    order.
    '''
    Z_2 = Integers(2) # Work mod 2
    order_dict = {}
    in_order = []
    place = 0 # Record the place in the ordering

    # Iterate over (Z/2)^2 in lexicographical order
    for a_0 in Z_2:
        for a_1 in Z_2:
            # Update the dict and list
            order_dict.update({(a_0, a_1) : place})
            in_order.append(vector([a_0, a_1]))
            place += 1

    return (order_dict, in_order)

def get_dim_single_skein(gamma):
    '''
    Takes an SL_2(Z) matrix gamma, and returns the dimension of the single skein
    part of the twisted torus, where gamma is an SL_2(Z)-matrix defining the
    twisting.

    The vector space in question is a quotient of C[X, Y]/(X^2 - 1 , Y^2 - 1) by
    some relations. The basis {1, X, Y, XY} is represented by vectors of length
    2 with Z/2-entries, i.e. X^aY^b is [a, b], and these are ordered
    lexicographically.

    The implementation is similar to get_relations: for all pairs of basis
    elements of C[X, Y]/(X^2 - 1, Y^2 - 1) we get the gamma-twisted commutators, then
    the corank of these relations is the required dimension.
    '''
    Z_2 = Integers(2) # Integers mod 2
    order_dict, in_order = order_lexi() # A dict and list to order the basis

    rels = [] # Ready to record the twisted commutator relations.

    # Walk thru pairs of basis elements (R, S) viewed as 1x2 Z/2 vectors (r, s)
    for r in in_order:
        for s in in_order:
            first_term = r + s # First term in the commutator is RS, i.e. r + s
            second_term = s*gamma.T + r # Second term is twisted by Gamma

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


def get_dim_estimates_empty(gamma, n, interactive_flag):
    '''
    Takes a matrix gamma and an integer n, and returns a list of estimates of
    the dimension of the skein module of the gamma-twisted torus, where the i-th
    element of the returned list is the estimated dimension at shell level i.

    If interactive_flag is true, prints verbosely to the command line.
    '''
    # Declare an indeterminate q.
    q = var('q')

    dimensions = [] #

    # Estimate the skein module dimension for each shell level.
    for shell_level in range(n+1):
        # For each shell level, compute #{lattice points}.
        N = (2*shell_level + 1)*(shell_level + 1) - shell_level
        if interactive_flag:
            print("Calculating relations for level %d (%d lattice points) ..." % (shell_level, N))
        relations = get_relations_empty(gamma, shell_level, order_lrtb)
        if interactive_flag:
            print("Found %d (non-independent) relations. Reducing ..." % len(relations))
        # Form a relation matrix, compute its pivots; the dimension estimate is the
        # co-rank.
        A = matrix(QQ['q'].fraction_field(), relations)
        pivots = A.pivots()
        dim_estimate = N - len(pivots)
        if interactive_flag:
            print("Dimension estimate for empty skein part at level %d: %d.\n\nVisualisation:\n" % (shell_level, dim_estimate))
            print_generators(shell_level, pivots, order_lrtb)

        dimensions.append(dim_estimate)

    return dimensions

def compute_and_write(M, shell_levels, path):
    '''
    A helper function for generate_raw_data, handles the subroutine of
    collating the results of dimension computations for M (up to shell_levels
    cutoff) and writes to the file at path, in append mode.
    '''
    # Get the dimension of the single skein, the esitmates for the empty skein
    # (not using interactive mode), and the sum of single dimension and the last
    # estimated empty dimension.
    dim_single = get_dim_single_skein(M)
    dim_estimates = get_dim_estimates_empty(M, shell_levels, False)
    dim_total = dim_single + dim_estimates[-1]

    # Write the relevant data to the output file.
    with open(path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([M.trace(), M[0, 0], M[0, 1], M[1, 0], M[1, 1], dim_single, dim_total] + dim_estimates)
    f.close()

def generate_raw_data(path, shell_levels):
    '''
    Generate several SL_2(Z) matrices, compute their skein dimension estimates,
    and write this data to a csv file.

    The matrices generated are grouped by absolute value of trace.
    '''

    #Trace 0 matrix
    S = matrix(ZZ, 2, [0, -1, 1, 0])

    #Trace 1 matrices
    M_0 = matrix(ZZ, 2, [1, -1, 1, 0])
    M_1 = matrix(ZZ, 2, [0, 1, -1, 1])

    low_trace = [S, M_0, M_1]

    #Generators for |trace| >= 2 matrices
    R = matrix(ZZ, 2, [1, 1, 0, 1])
    L = matrix(ZZ, 2, [1, 0, 1, 1])

    # Open a file for the raw data, and write a header.
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["trace", "a", "b", "c", "d", "single_dim", "total_dim"] + ["shell_{n}".format(n=i) for i in range(shell_levels + 1)])
    f.close()

    # Compute dimensions for the 3 matrices of low trace.
    for M in low_trace:
        compute_and_write(M, shell_levels, path)

    # Dimensions for the family of shears.
    for n in range(6):
        M = R**n
        compute_and_write(M, shell_levels, path)

    # Dimensions for words of length 2 in SL_2(Z)
    for m in range(1, 6):
        for n in range(6):
            M = (R**n)*(L**m)
            compute_and_write(M, shell_levels, path)

def write_dim_table(rawpath, outpath, shell_levels):
    '''
    Parses data produced by generate_raw_data, stored in the file rawpath, and
    gives this as a formatted table written to the file at outpath.
    '''
    df  = pd.read_csv(rawpath, dtype=str)

    #Calculate values that may be required for padding and alignment.
    max_tr_len = df["trace"].map(len).max()
    max_entry_len = df[["a", "b", "c", "d"]].applymap(len).values.max()
    max_empty_dim_len = df["shell_{n}".format(n=shell_levels)].map(len).max()
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
            print("{empty: >{max_empty_dim_len}s}\t".format(empty=row["shell_{n}".format(n=shell_levels)], max_empty_dim_len = max(max_empty_dim_len, len("EMPTY"))), end="", file=f)
            print("{tot: >{max_total_dim_len}s}\t\t\t".format(tot=row["total_dim"], max_total_dim_len=max(max_total_dim_len, len("TOTAL"))), end="", file=f)
            print("{ests}".format(ests=row[["shell_{n}".format(n=i) for i in range(shell_levels + 1)]].values.tolist()), file=f)
            print(" "*max(max_tr_len,  5) + "\t\t" + " ", end="", file=f)
            print("[{c: >{max_entry_len}s} {d: >{max_entry_len}s}]]".format(c=row["c"], d=row["d"], max_entry_len=max_entry_len), file=f)
            print("", file=f)
        f.close()

# Handle command-line flag
if len(sys.argv) < 2:
    print("Error, requires a mode choice (one of: i, g, w, gw) as a command-line argument")
    sys.exit(1)

choice = sys.argv[1]

# Interactive mode
if choice == "i":
    print("TWISTED TORUS SKEIN DIMENSION ESTIMATOR")
    print("Input an SL_2(Z)-matrix to define a twisted 3-torus.\nEnter a matix I, S, T, or enter the matrix\n\n[[a b]\n [c d]]\n\nas the string a b c d, and press return.")
    user_input = input().split(" ")

    # Allow users to specify the matrixes I, S, T just by letters
    if len(user_input) == 1:
        if user_input[0] == "I":
            print("Estimating dimensions for the matrix gamma = I...\n")
            gamma = matrix(ZZ, 2, [1, 0, 0, 1])
        elif user_input[0] == "S":
            print("Estimating dimensions for the matrix gamma = S...\n")
            gamma = matrix(ZZ, 2, [0, -1, 1, 0])
        elif user_input[0] == "T":
            print("Estimating dimensions for the matrix gamma = T...\n")
            gamma = matrix(ZZ, 2, [1, 1, 0, 1])
        elif user_input[0] == "-I":
            print("Estimating dimensions for the matrix gamma = -I...\n")
            gamma = matrix(ZZ, 2, [-1, 0, 0, -1])
        elif user_input[0] == "-S":
            print("Estimating dimensions for the matrix gamma = -S...\n")
            gamma = matrix(ZZ, 2, [0, 1, -1, 0])
        elif user_input[0] == "-T":
            print("Estimating dimensions for the matrix gamma = -T...\n")
            gamma = matrix(ZZ, 2, [-1, -1, 0, -1])
        # These are the only accepted special input matrices.
        else:
            print("Invalid input!")
            sys.exit(1)
    #Otherwise,  accept a list of 4 integers. Catch if not valid.
    elif len(user_input) != 4:
        print("Error! You did not enter 4 space-separated integers.")
        sys.exit(1)
    # Form a matrix from the user-input list, check it is a valid SL_2(Z) matrix
    else:
        # This matrix defines the twisted torus
        gamma = matrix(ZZ, 2, [int(i) for i in user_input])
        print("You have entered the matrix \n\n[[%s %s]\n [%s %s]]\n" % tuple(user_input))

        if gamma.determinant() != 1:
            print("Error: the data you entered is not an SL_2(Z) matrix, must have determinant 1.")
            sys.exit(1)

    # Having handled the user's choice of matrix, compute the dimension of the
    # single skein part.
    print("Dimension for single skein part: %d\n" % get_dim_single_skein(gamma))

    # Solicit the shell level to check up to from the user.
    n = int(input("Enter the number of shell levels to estimate dimension of empty skein part: "))

    # Get the dimension estimates, pass True interactive_flag to be verbose.
    get_dim_estimates_empty(gamma, n, True)

# Generation mode
elif choice == "g":
    path = "skein-dims-rawdata.csv"
    shell_levels = 6
    if len(sys.argv) >= 3:
        path = sys.argv[2]
        if len(sys.argv) >= 4:
            shell_levels = int(sys.argv[3])

    generate_raw_data(path, shell_levels)

# Presentation mode:
elif choice == "w":
    rawpath = "skein-dims-rawdata.csv"
    outpath = "skein-dims-printed.txt"
    shell_levels = 6
    if len(sys.argv) >= 3:
        rawpath = sys.argv[2]
        if len(sys.argv) >= 4:
            outpath = sys.argv[3]
            if len(sys.argv) >= 5:
                shell_levels = int(sys.argv[4])

    write_dim_table(rawpath, outpath, shell_levels)

# Generate-write mode:
elif choice == "gw":
    rawpath = "skein-dims-rawdata.csv"
    outpath = "skein-dims-printed.txt"
    shell_levels = 6
    if len(sys.argv) >= 3:
        rawpath = sys.argv[2]
        if len(sys.argv) >= 4:
            outpath = sys.argv[3]
            if len(sys.argv) >= 5:
                shell_levels = int(sys.argv[4])

    generate_raw_data(rawpath, shell_levels)
    write_dim_table(rawpath, outpath, shell_levels)

# Exit if an pinvalid mode choice is made.
else:
    print("Invalid choice of mode!")
    sys.exit(1)
