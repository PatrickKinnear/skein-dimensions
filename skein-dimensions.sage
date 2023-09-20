#!/usr/bin/env sage

'''
A  sxript to run the generate_raw_data function from skeinslib

USAGE Ensure SAGE_ROOT is stored in your PATH, and run

./skein-dimensions.sage

or run

sage skein-dimensions.sage

'''

import sys
from sage.all import *

load("skeinslib.sage")

generate_raw_data(max_seq_len=20)
