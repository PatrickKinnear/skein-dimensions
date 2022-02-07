# Skein Module Dimension Estimates

A project to estimate the dimensions of skein modules of twisted 3d tori.

Takes growing shells of lattice points representing generators, finds relations, and reduces.

If dimension estimates stabilise as the shells grow, this suggests a version of the diamond lemma holds and gives an upper bound on the dimension.

Written by Patrick Kinnear and Alisa Sheinkman.

## Usage

Ensure SAGE_ROOT is stored in your PATH, and run

./skein-dimensions.sage mode [rawpath [outpath]] [shell-level]

where mode is a string specifying the mode, one of:

(i) interactive
(g) generation
(w) write
(gw) generate-and-write

*Note* that w and gw modes require pandas. This can be installed by opening a
Sage shell and running

pip install pandas

*Requirements* Sage v9.2 or higher; Python v3.7 or higher.

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

## Description

A script to estimate the dimension of the skein module of the
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
