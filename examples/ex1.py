#!/bin/python

import pandas as pd
import numpy as np
import taguchi as tg

d = {
    'X1':[0,1,2],
    'X2':[1,2,3],
    'X3':[5,6,7],
    'X4':[1.,9.,22.],
    'X5':[-2,20,60],
    'X6':[-20,5,30],
}

n_levels = max([len(p) for p in d.values()])
n_params = len(d)

indexs = ['inputs','SN']
levels = range(0,n_levels)

midx = pd.MultiIndex.from_product([indexs,levels])

df = pd.DataFrame(np.zeros([len(midx),n_params]),columns=d.keys(),index=midx)
df.loc['inputs'] = pd.DataFrame(d).values

ta = tg.array(n_levels,n_params)

X = []
for n in range(0,n_params):
    X.append(df.loc['inputs'].iloc[ta[:,n],n])

np.array(X).T
