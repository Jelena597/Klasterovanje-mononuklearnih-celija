import pandas as pd
import numpy as np
from matplotlib import pyplot
import umap
from sklearn import cluster
from sklearn.decomposition import PCA

RANDOM_SEED = 42
DATA_LOCATION = 'data/interim/novi.csv'
COLORS = {
    'NK': 'cyan',
    'MC': 'yellow',
    'M+DC': 'blue',
    'T': 'green',
    'B': 'red',
}
DEFAULT_TYPES = {
    'NK': 0,
    'MC': 0,
    'M+DC': 0,
    'T': 0,
    'B': 0,
}
START_SAMPLE = 3000
SAMPLE_SIZE = 5000
NUMBER_OF_CLUSTERS = 8
EIGEN_TOL = 0.02
ALGORITHM = 'spectral'
# cluto
# cluto
# cluto
# cluto

np.random.seed(RANDOM_SEED)


def _read_data(data_location):
    return pd.read_csv(data_location, skiprows=START_SAMPLE, nrows=SAMPLE_SIZE)


def _reduce_dimensions_umap(data):
    umap_model = umap.UMAP(n_components=2, min_dist=0.5, metric='cosine', verbose=True,
                           random_state=RANDOM_SEED)
    umap_coords = umap_model.fit_transform(data)
    return umap_coords


def _cluster_data_spectral(data):
    return cluster.SpectralClustering(
        n_clusters=NUMBER_OF_CLUSTERS,
        random_state=RANDOM_SEED,
        affinity='nearest_neighbors',
        eigen_tol=EIGEN_TOL,
    ).fit(data)


def _cluster_data_kmeans(data):
    return cluster.KMeans(
        n_clusters=NUMBER_OF_CLUSTERS,
        random_state=RANDOM_SEED,
        verbose=True,
        algorithm='elkan',
    ).fit(data)


def _reduce_dimensions_pca(data, dim):
    return PCA(n_components=dim).fit(data).transform(data)


def _pie_chart(cluster_labels, original_labels):
    cluster_predictions = {i: DEFAULT_TYPES.copy() for i in range(NUMBER_OF_CLUSTERS)}
    for cluster_label, original_label in zip(cluster_labels, original_labels):
        cluster_predictions[cluster_label][original_label] += 1
    print(cluster_predictions)
    ncols = 3
    fig1, ax1 = pyplot.subplots(nrows=3, ncols=ncols)
    for i in range(NUMBER_OF_CLUSTERS):
        ax1[int(i/ncols), i%ncols].pie(
            cluster_predictions[i].values(),
            # labels=cluster_predictions[i].keys(),
            autopct='%1.1f%%',
            startangle=90,
            colors=COLORS.values(),
        )
        ax1[int(i/ncols), i%ncols].axis('equal')
        ax1[int(i/ncols), i%ncols].set_title(f'Broj celija: {sum(cluster_predictions[i].values())}')
    pyplot.legend(COLORS.keys())
    pyplot.savefig(f'results/pie_{ALGORITHM}_clusters_{NUMBER_OF_CLUSTERS}_size_{SAMPLE_SIZE}_eigen_{EIGEN_TOL}.png')
    pyplot.show()


def main():
    data = _read_data(DATA_LOCATION)
    data_table_values = data.iloc[:, 3:]

    # Samo umap
    umap_coords = _reduce_dimensions_umap(data_table_values)
    fig, ax = pyplot.subplots(2)
    ax[0].scatter(umap_coords[:, 0], umap_coords[:, 1], s=10, alpha=.5, marker="X",
                  c=data.iloc[:, 2].map(COLORS))

    data_pca = _reduce_dimensions_pca(data_table_values, 50)

#     # umap_pca_coords = _reduce_dimensions_umap(data_pca)
#     # ax[1].scatter(umap_pca_coords[:, 0], umap_pca_coords[:, 1], s=10, alpha=.5, marker="X",
#     #               c=data.iloc[:, 2].map(COLORS))
#
#     if ALGORITHM == 'spectral':
#         clustering = _cluster_data_spectral(data_table_values)
#         ax[1].scatter(umap_coords[:, 0], umap_coords[:, 1], s=10, alpha=.5, marker="X",
#                       c=clustering.labels_)
#         # _pie_chart(clustering.labels_, data.iloc[:, 2])
#         pyplot.savefig(f'results/{ALGORITHM}_clusters_{NUMBER_OF_CLUSTERS}_eigen_{EIGEN_TOL}.png')

    if ALGORITHM == 'kmeans':
        clustering = _cluster_data_kmeans(data_table_values)
        ax[1].scatter(umap_coords[:, 0], umap_coords[:, 1], s=10, alpha=.5, marker="X",
                      c=clustering.labels_)
        # _pie_chart(clustering.labels_, data.iloc[:, 2])
        pyplot.savefig(f'results/{ALGORITHM}_clusters_{NUMBER_OF_CLUSTERS}')

    pyplot.show()


if __name__ == '__main__':
    main()
