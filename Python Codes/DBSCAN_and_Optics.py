import time
import warnings
from itertools import cycle, islice

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn import cluster, datasets, mixture
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler

## Check for NaN values in the input file.
file_name = "Cluster_DIFFUSION.txt"
data_file = np.loadtxt(file_name)
np.random.shuffle(data_file)
print('Total data in file :', len(data_file), data_file.shape)
print('cluster density and average deviation angle data file:')
print(data_file[:10])
df = pd.DataFrame(data_file)
print('df.shape: ', df.shape)
eighty_pct = 0.8*df.shape[0]


X = df.loc[:, 0:1] # 0th and 1st columns
y = df.loc[:, 2] # 2nd columnX = 

### Scale data (two possibilities - choose one)
# scalar = MinMaxScaler()
# X = scalar.fit_transform(X)
scalar = StandardScaler()
X = scalar.fit_transform(X)
print('Scaled:', X[:8])
blobs = X
varied = X

# Plot input PoCA points
# colors = [plt.cm.jet(each) for each in np.linspace(0, 1, 5)]
color = 'green'
plt.figure(figsize=(9, 6))
plot_sample = 10000
# Following plotting instruction does not work in Jupyter notebook
for xp, yp, zp in zip(X[:, 0], X[:, 1], y):
    plt.plot(xp, yp, c = color, marker='o', markerfacecolor='none', markersize='4')
# Following works, instead.
# xp = X_train.iloc[:, 0]
# yp = X_train.iloc[:, 1]
# plt.plot(xp, yp, 'o', markersize=4)
plt.grid()
plt.xlabel('Average cluster density (per $25.0 mm^2$)', fontsize=17)
plt.ylabel('Average deviation angle (deg)', fontsize=17)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title('Deviation angle vs cluster density')
plt.savefig('ClusterAlgos_ScatterPlot_CompareClusterAlgos.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()



# ============
# Set up cluster parameters
# ============
plt.figure(figsize=(9 * 2 + 3, 13))
plt.subplots_adjust(
    left=0.02, right=0.98, bottom=0.001, top=0.95, wspace=0.05, hspace=0.01
)

plot_num = 1

default_base = {
    "quantile": 0.3,
    "eps": 0.3,
    "damping": 0.9,
    "preference": -200,
    "n_neighbors": 3,
    "n_clusters": 5,
    "n_init": 'auto',
    "min_samples": 7,
    "xi": 0.05,
    "min_cluster_size": 0.1,
    "allow_single_cluster": True,
    "random_state": 42,
}

# update parameters with dataset-specific values
params = default_base.copy()
#params.update(algo_params)

# estimate bandwidth for mean shift
bandwidth = cluster.estimate_bandwidth(X, quantile=params["quantile"])

# connectivity matrix for structured Ward
connectivity = kneighbors_graph(
    X, n_neighbors=params["n_neighbors"], include_self=False
)
# make connectivity symmetric
connectivity = 0.5 * (connectivity + connectivity.T)

# ============
# Create cluster objects
# ============

dbscan = cluster.DBSCAN(eps=params["eps"])
optics = cluster.OPTICS(
    min_samples=params["min_samples"],
    xi=params["xi"],
    min_cluster_size=params["min_cluster_size"],
)

clustering_algorithms = (
    ("DBSCAN", dbscan),
    ("OPTICS", optics),
    
)

for name, algorithm in clustering_algorithms:
    t0 = time.time()

    # catch warnings related to kneighbors_graph
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="the number of connected components of the "
            + "connectivity matrix is [0-9]{1,2}"
            + " > 1. Completing it to avoid stopping the tree early.",
            category=UserWarning,
        )
        warnings.filterwarnings(
            "ignore",
            message="Graph is not fully connected, spectral embedding"
            + " may not work as expected.",
            category=UserWarning,
        )
        algorithm.fit(X)

    t1 = time.time()
    if hasattr(algorithm, "labels_"):
       y_pred = algorithm.labels_.astype(int)
    else:
       y_pred = algorithm.predict(X)

    plt.subplot(1, len(clustering_algorithms), plot_num)

    plt.title(name, size=18)
    colors = np.array(
        list(
            islice(
                cycle(
                    [
                        "#377eb8",
                        "#ff7f00",
                        "#4daf4a",
                        "#f781bf",
                        "#a65628",
                        "#984ea3",
                        "#999999",
                        "#e41a1c",
                        "#dede00",
                    ]
                ),
                int(max(y_pred) + 1),
            )
        )
    )
    colors = np.append(colors, ["#000000"])
    plt.scatter(X[:, 0], X[:, 1], s=10, color=colors[y_pred])

    plt.xlim(-2.5, 2.5)
    plt.ylim(-2.5, 2.5)
    plt.xticks(())
    plt.yticks(())
    plt.text(
        0.99,
        0.01,
        ("%.2fs" % (t1 - t0)).lstrip("0"),
        transform=plt.gca().transAxes,
        size=15,
        horizontalalignment="right",
    )
    plot_num += 1

plt.savefig('DBSCAN_and_Optics.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()



