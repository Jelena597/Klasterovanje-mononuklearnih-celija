import pandas as pd
import numpy as np
from sklearn.cluster import KMeans, SpectralClustering, DBSCAN
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.preprocessing import normalize
import anndata
import seaborn as sb
from sklearn.neighbors import kneighbors_graph, DistanceMetric
from sklearn.metrics import silhouette_score


def add_neighbors_and_umapCoords(adata, n_neighs, metric_name, min_d):
    ad1 = sc.AnnData(X=adata.X, var=adata.var, obs=adata.obs)

    sc.pp.neighbors(ad1, n_neighbors=n_neighs, use_rep="X", metric=metric_name)
    sc.tl.umap(ad1, min_dist=min_d, n_components=2)

    return ad1, metric_name

def cluster(ad, algorithm, n_cls, metricName):
    labels = []
    if algorithm.lower() == 'kmeans':
        clustering = KMeans(n_clusters=n_cls).fit(ad.X)
        labels = clustering.labels_
    elif algorithm.lower() == 'spectral':
        ds = ad.obsp["distances"]
        ds = pd.DataFrame.sparse.from_spmatrix(ds)
        ds_s = (ds + ds.T)/2
        clustering = SpectralClustering(n_clusters=n_cls, affinity='precomputed').fit(ds_s)
        labels = clustering.labels_
    y = pd.Series(labels)
    labeling(y, ad, metricName)

def labeling(y, ad, metricName):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,4))
    fig.suptitle(metricName)

    colors1 = []
    colors2 = []
    exacts = []
    expected = []

    for index, value in y.items():
        if value == 0:
            colors1.append('cl1')
        elif value == 1:
            colors1.append('cl2')
        elif value == 2:
            colors1.append('cl3')
        elif value == 3:
            colors1.append('cl4')
        elif value == 4:
            colors1.append('cl5')
        elif value == 5:
            colors1.append('cl6')
        elif value == 6:
            colors1.append('cl7')
        elif value == 7:
            colors1.append('cl8')
        elif value == 8:
            colors1.append('cl9')
        elif value == 9:
            colors1.append('cl10')
        elif value == 10:
            colors1.append('cl11')
        elif value == 11:
            colors1.append('cl12')
        elif value == 12:
            colors1.append('cl13')
        elif value == 13:
            colors1.append('cl14')
        elif value == 14:
            colors1.append('cl15')


    for i in range(1786):
        exacts.append('B')

    for i in range(2445):
        exacts.append('MDC')

    for i in range(856):
        exacts.append('MC')

    for i in range(309):
        exacts.append('NK')

    for i in range(325):
        exacts.append('iNKT')

    for i in range(382):
        exacts.append('MaiT')

    for i in range(4486):
        exacts.append('P5T')

    for i in range(1247):
        exacts.append('T4_1')

    for i in range(222):
        exacts.append('T4_2')

    for i in range(965):
        exacts.append('T4_3')

    for i in range(435):
        exacts.append('T4_4')

    for i in range(310):
        exacts.append('T8')

    for i in range(4753):
        exacts.append('T8_metha')

    for i in range(284):
        exacts.append('Vd1')

    for i in range(204):
        exacts.append('Vd2')

    exacts = pd.Series(exacts)

    names = pd.Series(ad1.obs_names)
    names_connect = pd.concat((names, y, exacts),axis=1)
    names_connect = names_connect.rename({0:'Cells', 1:'Cluster', 2:'Names'},axis=1)
    sort_cells = names_connect.sort_values(['Names'])

    for i, v in sort_cells.iterrows():
        if v['Names'] == 'B':
            colors2.append('B')
        elif v['Names'] == 'MDC':
            colors2.append('MDC')
        elif v['Names'] == 'MC':
            colors2.append('MC')
        elif v['Names'] == 'NK':
            colors2.append('NK')
        elif v['Names'] == 'iNKT':
            colors2.append('T1')
        elif v['Names'] == 'MaiT':
            colors2.append('T23')
        elif v['Names'] == 'P5T':
            colors2.append('T23')
        elif v['Names'] == 'T4_1':
            colors2.append('T4')
        elif v['Names'] == 'T4_2':
            colors2.append('T5678')
        elif v['Names'] == 'T4_3':
            colors2.append('T5678')
        elif v['Names'] == 'T4_4':
            colors2.append('T5678')
        elif v['Names'] == 'T8':
            colors2.append('T5678')
        elif v['Names'] == 'T8_metha':
            colors2.append('T9_10')
        elif v['Names'] == 'Vd1':
            colors2.append('T9_10')
        elif v['Names'] == 'Vd2':
            colors2.append('T11')

    s = silhouette_score(ad.X, y)
    print(s)

    ad.obs['ClassFound'] = colors1
    found = sc.pl.umap(ad, color="ClassFound", ax = ax1, show=False)

    ad.obs['ClassReal'] = colors2

    real = sc.pl.umap(ad, color="ClassReal", ax = ax2, show=False)


#specify the path
adata = sc.read_csv('C:\\Users\\User\Desktop\seminarski.ip2\GEO\\klaster\\Sve.csv')

#example
ad1, metric = add_neighbors_and_umapCoords(adata, 110, 'cosine', 0.1)
cluster(ad1, 'kmeans', 6, metric)

