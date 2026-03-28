import numpy as np
import pyemma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pyemma.util.contexts import settings
from sklearn.preprocessing import StandardScaler
from scipy.ndimage import gaussian_filter
from scipy.ndimage import minimum_filter
from pyemma.msm import tpt

def main():

    # =========================
    # PARAMETERS
    # =========================

    feature_file = "tica_msm_features_signif.dat"

    tica_lag = 100
    msm_lag = 25

    n_clusters = 100
    n_metastable = 4

    # =========================
    # LOAD FEATURES
    # =========================

    print("Loading features...")

    data = np.loadtxt(feature_file)

    X = data[2420:,1:]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    print("Feature matrix:", X_scaled.shape)

    # =========================
    # tICA
    # =========================

    print("Running tICA...")

    tica = pyemma.coordinates.tica(
        X_scaled,
        lag=tica_lag,
        dim=5,
        kinetic_map=True
    )

    Y = tica.get_output()[0]

    print("tICA output shape:", Y.shape)

    # =========================
    # CLUSTERING
    # =========================

    print("Clustering...")

    cluster = pyemma.coordinates.cluster_kmeans(
        Y,
        k=n_clusters,
        max_iter=1000
    )

    dtrajs = cluster.dtrajs

    # =========================
    # IMPLIED TIMESCALES
    # =========================

    print("Computing implied timescales...")
    its = pyemma.msm.its(dtrajs, lags=[2,4,6,8,10,12,14,16,18,20,25,30,40,45,50,100,150], nits=17, errors='bayes')

    pyemma.plots.plot_implied_timescales(its)
    plt.xlim(2,50)
    plt.xlabel('lagtime/steps',fontsize=24,fontweight='bold')
    plt.ylabel('timescale/steps',fontsize=24,fontweight='bold')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # Axis styling
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    plt.tight_layout()
    plt.savefig('sim1-it-lag-filter.png')
    plt.show()

if __name__ == "__main__":
    main()

