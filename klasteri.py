import pandas as pd
import numpy as np
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics import v_measure_score
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import seaborn as sb
from scipy.sparse import csgraph

adata = sc.read_csv('../data/interim/geni_shuffled.csv')

ad1 = sc.AnnData(X=adata.X, var=adata.var, obs=adata.obs)

sc.pp.neighbors(ad1, n_neighbors=200, use_rep="X", metric='minkowski')
#print(ad1.obsp)
ds = ad1.obsp["distances"]
ds = pd.DataFrame.sparse.from_spmatrix(ds)

ds_s = (ds + ds.T)/2

#print(ds_s)

sc.tl.umap(ad1, min_dist=0.1, n_components=2)

clustering = SpectralClustering(n_clusters=4, affinity='precomputed').fit(ds_s)
labels = clustering.labels_
y = pd.Series(labels)

colors = []

for index, value in y.items():
    if value == 0:
        colors.append(0)
    elif value == 1:
        colors.append(1)
    elif value == 2:
        colors.append(2)
    else:
        colors.append(3)

label = []

for i in range(1786):
    label.append('B')

for i in range(2445):
    label.append('MDC')

for i in range(856):
    label.append('MC')

for i in range(309):
    label.append('NK')

label = pd.Series(label)


B_0 = 0
B_1 = 0
B_2 = 0
B_3 = 0

MDC_0 = 0
MDC_1 = 0
MDC_2 = 0
MDC_3 = 0

MC_0 = 0
MC_1 = 0
MC_2 = 0
MC_3 = 0

NK_0 = 0
NK_1 = 0
NK_2 = 0
NK_3 = 0


names = pd.Series(ad1.obs_names)
names_connect = pd.concat((names, y, label),axis=1)
names_connect = names_connect.rename({0:'Cells', 1:'Cluster', 2:'Names'},axis=1)
sort_cells = names_connect.sort_values(['Names'])

for i, v in sort_cells.iterrows():
    if v['Names'] == 'B':
        if v['Cluster'] == 0:
            B_0 += 1
        elif v['Cluster'] == 1:
            B_1 += 1
        elif v['Cluster'] == 2:
            B_2 += 1
        else:
            B_3 += 1
    elif v['Names'] == 'MDC':
        if v['Cluster'] == 0:
            MDC_0 += 1
        elif v['Cluster'] == 1:
            MDC_1 += 1
        elif v['Cluster'] == 2:
            MDC_2 += 1
        else:
            MDC_3 += 1
    elif v['Names'] == 'MC':
        if v['Cluster'] == 0:
            MC_0 += 1
        elif v['Cluster'] == 1:
            MC_1 += 1
        elif v['Cluster'] == 2:
            MC_2 += 1
        else:
            MC_3 += 1
    else:
        if v['Cluster'] == 0:
            NK_0 += 1
        elif v['Cluster'] == 1:
            NK_1 += 1
        elif v['Cluster'] == 2:
            NK_2 += 1
        else:
            NK_3 += 1


ad1.obs['Class'] = colors
sc.pl.umap(ad1, color="Class")

