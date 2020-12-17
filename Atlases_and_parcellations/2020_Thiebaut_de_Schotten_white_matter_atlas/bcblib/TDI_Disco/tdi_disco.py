""" This script aims to estimate the optimal number of healthy subject
tractographies we need to have an optimal result with Disconnectome Maps.
"""
import sys

import pandas as pd
import numpy as np

groups = pd.read_csv(sys.argv[1], header=None, index_col=None)

print(str(np.array(groups.shape)))

for c in groups.columns.values:
    print("line: " + str(c))
    for l in groups.index.values:
        print(groups[l][c])
