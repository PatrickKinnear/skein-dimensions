#!/usr/bin/env sage

'''
USAGE Ensure SAGE_ROOT is stored in your PATH, and run

./skein-dimensions.sage mode [rawpath [outpath]] [shell-level]

where mode is a string specifying the mode, one of:

(i) interactive
(g) generation
(w) write
(gw) generate-and-write

*Requirements* Sage v9.2 or higher; Python v3.7 or higher.

*Note* that w and gw modes require pandas. This can be installed by opening a
Sage shell and running

pip install pandas

or by running

sage -pip install pandas

at the command line.

Most of the functionality is implemented in the library skeinslib.sage, which
should be in the same directory as this script.

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
skein-dims-rawdata.csv. Can also specify shell level (defaults to 11).

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

import sys
from sage.all import *

load("skeinslib.sage")

# Handle command-line flag
if len(sys.argv) < 2:
    print("Error, requires a mode choice (one of: i, g, w, gw) as a command-line argument")
    sys.exit(1)

choice = sys.argv[1]

shell_levels = 11

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
    compute_reduced_matrix(gamma, n, True)

# Generation mode
elif choice == "g":
    path = "skein-dims-rawdata.csv"
    if len(sys.argv) >= 3:
        path = sys.argv[2]
        if len(sys.argv) >= 4:
            shell_levels = int(sys.argv[3])

    generate_raw_data(path, shell_levels)

# Presentation mode:
elif choice == "w":
    import pandas as pd
    rawpath = "skein-dims-rawdata.csv"
    outpath = "skein-dims-printed.txt"
    if len(sys.argv) >= 3:
        rawpath = sys.argv[2]
        if len(sys.argv) >= 4:
            outpath = sys.argv[3]
            if len(sys.argv) >= 5:
                shell_levels = int(sys.argv[4])

    write_dim_table(rawpath, outpath, shell_levels)

# Generate-write mode:
elif choice == "gw":
    import pandas as pd
    rawpath = "skein-dims-rawdata.csv"
    outpath = "skein-dims-printed.txt"
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
