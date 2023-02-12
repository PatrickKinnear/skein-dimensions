'''
DESCRIPTION Some  functionality for tabulating dimensions of skein modules of
mapping tori if T^2. The formulae are given in a paper of P. Kinnear, and this
code simply automates the computations.

USAGE In a sage interactive session or a sage script, load these functions using

load("skeinslib.sage")

*Requirements* Sage v9.2 or higher; Python v3.7 or higher.

*Note* that the function write_dim_table required pandas imported as pd. This
can be installed by opening a Sage shell and running

pip install pandas

or by running

sage -pip install pandas

at the command line.
'''

import csv
import itertools
import os
import pandas as pd

from sage.all import *

def get_dim_single_skein(gamma):
    '''
    Takes an SL_2(Z) matrix gamma, and returns the dimension of the single skein
    part of the twisted torus, where gamma is an SL_2(Z)-matrix defining the
    twisting.

    The dimension depends only on the parity of gamma, and this code simply
    checks this. See the accompanying paper for the derivation of the dims.
    '''
    I = matrix(ZZ, 2, [1, 0, 0, 1])
    E_plus_2 = matrix(ZZ, 2, [1, 1, 1, 0]) # The matrix E_plus mod 2
    E_minus_2 = matrix(ZZ, 2, [0, 1, 1, 1]) # The matrix E_minus mod 2

    if gamma % 2 == I:
        return 4
    elif gamma % 2 == E_plus_2 or gamma % 2 == E_minus_2:
        return 1
    else:
        return 2


def weyl_action(x):
    '''
    Implements taking coinvariants for the negation action of Z/2Z on quotients
    Z/xZ. Given x defining Z/xZ, return the number of elements of the space of
    coinvariants.
    '''
    if x % 2 == 0:
        return (x + 2)/2
    else:
        return (x + 1)/2

def get_dim_empty_skein(gamma):
    '''
    Takes an SL_2(Z) matrix gamma, and returns the dimensions of the empty skein
    part of the twisted torus, where gamma is an SL_2(Z)-matrix defining the
    twisting.

    The dimension depends only the invariant factors of I +/- gamma, acting on
    the lattice Z^2, and this is obtained by smith normal form.

    The empty skein part is a direct sum of two components, correpsonding to
    I - gamma and I + gamma. A list of the dimension of each summand is returned
    '''
    I = matrix(ZZ, 2, [1, 0, 0, 1])

    D_minus, U_minus, V_minus = (I - gamma).smith_form()
    a_minus = D_minus.diagonal()

    D_plus, U_plus, V_plus = (I + gamma).smith_form()
    a_plus = D_plus.diagonal()

    n_minus = [weyl_action(a) for a in a_minus]
    n_plus = [weyl_action(a) for a in a_plus]

    return [prod(n_minus), prod(n_plus)]

def compute_and_write(sequence, gamma, output_path, cache_path):
    '''
    A helper function for generate_raw_data, handles the subroutine of
    collating the results of dimension computations for M and writes to the file
    at output_path. Also maintains a persistent cache.
    '''
    # Get the dimension of the single skein.
    dim_single = get_dim_single_skein(gamma)

    # Dimension of empty skein part.
    dim_empty_minus, dim_empty_plus = get_dim_empty_skein(gamma)

    # Total
    dim_total = dim_single + dim_empty_minus + dim_empty_plus

    # Write the relevant data to the output files.
    with open(output_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([gamma.trace(), gamma[0, 0], gamma[0, 1], gamma[1, 0], gamma[1, 1], dim_single, dim_empty_minus, dim_empty_plus, dim_total] + list(sequence))
    f.close()

    # Write the computed sequence to the persistent cache file.
    with open(cache_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(list(sequence))

    return None

def compute_write_low_trace(output_path, cache_path):
    '''
    Perform the computations for matrices of trace with abs val < 2, and write
    to the file at path.
    '''
    #Trace 0 matrix
    S = matrix(ZZ, 2, [0, -1, 1, 0])

    #Trace 1 matrices
    E_plus = matrix(ZZ, 2, [1, -1, 1, 0])
    E_minus = matrix(ZZ, 2, [0, 1, -1, 1])

    low_trace = [S, E_plus, E_minus]

    # Compute dimensions for the 3 matrices of low trace.
    for M in low_trace:
        compute_and_write([], M, output_path=output_path, cache_path=cache_path)

    return None

def compute_write_from_seq(sequence, output_path, cache_path):
    '''
    Performs the dimension computations for an SL_2(Z) matrix of form
        R^{a_1}L^{a_2}...(R or L)^{a_k}
    given a sequence (a_k), and writes to path.

    Here R = [[1, 1], [0, 1]], L = [[1, 0], [1, 1]] so that R^n is
    [[1, n], [0, 1]] and this is how we implement the exponentiation (similar
    for L).
    '''
    M = matrix(ZZ, 2, [1, 0, 0, 1])

    # Build up the word indexed by this sequence
    for i in range(len(sequence)):
        if i % 2 == 0:
            #Multiply by R^i
            M = M*matrix(ZZ, 2, [1, sequence[i], 0, 1])
        else:
            #Multiply by L^i
            M = M*matrix(ZZ, 2, [1, 0, sequence[i], 1])

    compute_and_write(sequence=sequence, gamma=M, output_path=output_path, cache_path=cache_path)

    return None

def seq_has_been_checked(seq, cache):
    '''
    Checks if the given sequence, or any rotation thereof, is in the cache of
    already-checked sequences.
    '''
    for i in range(len(seq)):
        seq = seq[1:] + seq[:1] # Rotate the sequence once.
        if seq in cache: # Check cache membership.
            return True

    return False

def generate_raw_data(append=False, output_path="skeindims-rawdata.csv", cache_path="seq-cache.csv", max_seq_len=6, nshears=11):
    '''
    Generate several SL_2(Z) matrices, compute their skein dimensions, and write
    this data to a csv file.

    Maintains a cache of previously checked sequences. In append mode, this is
    loaded from a csv and the new dimension data is appended to an existing
    output file.
    '''
    half_max_seq_len = max_seq_len//2

    if not append:
        # Open a file for the raw data, and write a header.
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["trace", "a", "b", "c", "d", "single_dim", "dim_HH_0_-", "dim_HH_0_+", "total_dim"] + ["seq_{n}".format(n=i) for i in range(2*half_max_seq_len)])
        # Reset the cache
        if os.path.exists(cache_path):
            os.remove(cache_path)

    # Dimensions for low trace matrices.
    compute_write_low_trace(output_path=output_path, cache_path=cache_path)

    # Dimensions for the family of shears (|trace = 2|).
    for n in range(nshears):
        compute_write_from_seq(sequence=[n], output_path=output_path, cache_path=cache_path)

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
                compute_write_from_seq(sequence=sequence, output_path=output_path, cache_path=cache_path)
                cache.append(sequence)

    return None

def write_dim_table(rawpath="skeindims-rawdata.csv", outpath="skeindims-formatted.txt"):
    '''
    Parses data produced by generate_raw_data, stored in the file rawpath, and
    gives this as a formatted table written to the file at outpath.
    '''
    df  = pd.read_csv(rawpath, dtype=str)

    #Calculate values that may be required for padding and alignment.
    max_tr_len = df["trace"].map(len).max()
    max_entry_len = df[["a", "b", "c", "d"]].applymap(len).values.max()
    max_HH_0_min = df["dim_HH_0_-"].map(len).max()
    max_HH_0_pls = df["dim_HH_0_+"].map(len).max()
    max_total_dim_len = df["total_dim"].map(len).max()

    mat_line_length = 5 + max_entry_len*2

    with open(outpath, "w") as f:
        #Table header
        print("TRACE\t\tMATRIX" + " "*(mat_line_length-len("MATRIX")) + "\t\tSINGLE\tHH_0_-\tHH_0_+\tTOTAL", file=f)

        #Iterate through the entries and print rows.
        for idx, row in df.iterrows():
            print("{tr: >{max_tr_len}s}\t\t".format(tr=row["trace"], max_tr_len=max(max_tr_len, len("TRACE"))), end="", file=f)
            print("[[{a: >{max_entry_len}s} {b: >{max_entry_len}s}] \t\t".format(a=row["a"], b=row["b"], max_entry_len=max_entry_len), end="", file=f)
            print("{sing: >6s}\t".format(sing=row["single_dim"]), end="", file=f)
            print("{HH_0_min: >6s}\t".format(HH_0_min=row["dim_HH_0_-"]), end="", file=f)
            print("{HH_0_pls: >6s}\t".format(HH_0_pls=row["dim_HH_0_+"]), end="", file=f)
            print("{tot: >{max_total_dim_len}s}\t\t\t".format(tot=row["total_dim"], max_total_dim_len=max(max_total_dim_len, len("TOTAL"))), end="\n", file=f)
            print(" "*max(max_tr_len,  5) + "\t\t" + " ", end="", file=f)
            print("[{c: >{max_entry_len}s} {d: >{max_entry_len}s}]]".format(c=row["c"], d=row["d"], max_entry_len=max_entry_len), file=f)
            print("", file=f)
        f.close()

    return None
