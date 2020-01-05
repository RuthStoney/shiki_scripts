#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 14:01:42 2020

@author: ruthstoney
"""

import pandas as pd, seaborn as sns
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from sklearn.datasets import load_iris
sns.set(font="monospace")

iris = load_iris()
X, y = iris.data, iris.target
Df = pd.DataFrame(X, index = ["iris_%d" % (i) for i in range(X.shape[0])], columns = iris.feature_names)

DF_corr = Df.T.corr()
DF_dism = 1 - DF_corr   # distance matrix
DF_dism = dist2.copy()
linkage = hc.linkage(sp.distance.squareform(DF_dism), method='average')
sns.clustermap(DF_dism, row_linkage=linkage, col_linkage=linkage)
