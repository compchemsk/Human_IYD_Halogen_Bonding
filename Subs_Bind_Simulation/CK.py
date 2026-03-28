import numpy as np
import pyemma
import pyemma.plots as mplt
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

def main():

    # =========================
    # PARAMETERS
    # =========================

    feature_file = "tica_msm_features_signif.dat"

    lag = 30
    n_clusters = 6
    n_metastable = 6

    # =========================
    # LOAD FEATURES
    # =========================

    print("Loading features...")

    data = np.loadtxt(feature_file)

    # remove equilibration frames
    X = data[2420:,1:]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    print("Feature matrix shape:", X_scaled.shape)

    # =========================
    # tICA
    # =========================

    print("Running tICA...")

    tica = pyemma.coordinates.tica(
        X_scaled,
        lag=lag,
        dim=5,
        kinetic_map=True
    )

    Y = tica.get_output()[0]

    print("tICA output shape:", Y.shape)

    # =========================
    # CLUSTERING
    # =========================

    print("Running k-means clustering...")

    cluster = pyemma.coordinates.cluster_kmeans(
        Y,
        k=n_clusters,
        max_iter=10000
    )

    print("Number of clusters:", cluster.n_clusters)

    # =========================
    # BUILD MSM
    # =========================

    print("Estimating MSM...")

    msm = pyemma.msm.estimate_markov_model(
        cluster.dtrajs,
        lag=lag
    )

    print("MSM built.")
    print("Number of MSM states:", msm.nstates)

    # =========================
    # CK TEST
    # =========================

    print("Running Chapman-Kolmogorov test...")

    ck = msm.cktest(6)
    print(ck.predictions)
    print(ck.estimates)

    # =========================
    # CK PLOT
    # =========================
    import matplotlib.figure

# Patch legend compatibility for PyEMMA
    _old_legend = matplotlib.figure.Figure.legend

    def _new_legend(self, *args, **kwargs):
        if len(args) == 3:
            handles, labels, loc = args
            return _old_legend(self, handles, labels, loc=loc, **kwargs)
        return _old_legend(self, *args, **kwargs)

    matplotlib.figure.Figure.legend = _new_legend


    import pyemma.plots as mplt
    import matplotlib.pyplot as plt

    fig, axes = mplt.plot_cktest(ck)
# Set x-axis limit from 0 to 50
    for ax_row in axes:
        for ax in ax_row:
            ax.set_xlim(0, 50)

    plt.tight_layout()
    plt.savefig('CK-test-sim1-it.png',dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
