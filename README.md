# Skein Module Dimensions

A project to calculate the dimensions of skein modules of mapping tori of T^2, written by Patrick Kinnear.

The formulae are given in a paper of P. Kinnear, and this code simply automates the computations. Most code is found in the file `skeinslib.sage`. The main functions of interest are

`get_dim_single_skein`

and

`get_dim_empty_skein`

which take an SL_2(Z)-matrix `gamma` which defines a mapping torus, and return the dimension of certain direct summands of its Kauffman bracket skein module (for generic q). The total dimension is the sum of the values returned by these functions.

There are additional functions for helping to tabulate, store and display the calculations for conjugacy classes in SL_2(Z). These computations are stored in the `data` folder. The initial data in this folder, for 132636 conjugacy classes, was generated on 4th November 2022 by P. Kinnear; further updates will be logged by date they are committed to this repo.

*Note* that the code for writing a nicely-formatted table requires pandas. This can be installed by opening a
Sage shell and running

`pip install pandas`

or by running

`sage -pip install pandas`

at the command line.

*Requirements* Sage v9.2 or higher; Python v3.7 or higher.
