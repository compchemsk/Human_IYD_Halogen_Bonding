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

    lag = 30

    n_clusters = 1000
    n_metastable = 8

    # =========================
    # FEATURE HEADER
    # =========================

    header = ['frame', 'RMSD_helix_A', 'RMSD_MIA', 'd_MIA_FRA', 'I181_halo_binary', 'MIA_Alpha5B', 'd_MIA_FRB', 'RMSD_helix_B', 
            'waters_MIA_5A', 'MIB_Alpha5A', '104_hydro_binary', 'RMSD_MIB', 'MIA_Alpha5A', '244*_hydro_binary', 
            '248*_hydro_binary', 'loopA_alpha5B_dist', '240*_hydro_binary', '242*_hydro_binary', 'MIB_Alpha5B', 
            'waters_helix_5A', 'Y184_halo_binary', 'waters_MIB_5A', 'loopB_alpha5A_dist', '177_hydro_binary', 
            'T171_halo_binary', '243*_hydro_binary', 'waters_helix_5B', '258*_hydro_binary', 'R248*_halo_binary'] 
    
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
        lag=lag,
        dim=5,
        kinetic_map=True
    )

    Y = tica.get_output()[0]

    print("tICA output shape:", Y.shape)

    # =========================
    # tICA TIMESERIES PLOT
    # =========================

    plt.figure(figsize=(10,5))

    plt.subplot(2,1,1)
    plt.plot(Y[:,0])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("tIC1",fontsize=18,fontweight='bold')
    plt.xlim(0,2000)
    plt.ylim(-2.5,2.5)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)

    plt.subplot(2,1,2)
    plt.plot(Y[:,1])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("tIC2",fontsize=18,fontweight='bold')
    plt.xlabel("Time (1 ns)",fontsize=18,fontweight='bold')
    plt.xlim(0,2000)
    plt.ylim(-2.5,2.5)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)

    plt.tight_layout()
    plt.savefig('tic1_tic2_sim1_it.png',dpi=300)   
    plt.show()


    # =========================
    # FREE ENERGY SURFACE
    # =========================

    print("Computing Free Energy Surface...")

    z,x,y = np.histogram2d(Y[:,0],Y[:,1],bins=50)

    P = z / np.sum(z)

    P[P==0] = 1e-12

    F = -np.log(P)

    F = F - np.min(F)

    F = gaussian_filter(F, sigma=1.0)

    np.savetxt('fes_tic1_tic2_sim1_it.dat',Y)
    extent = [-2.5,2.5,-2.5,2.5]

    plt.figure(figsize=(8,6))

    contour = plt.contourf(
        F.T,
        levels=30,
        cmap=cm.jet,
        extent=extent,
        origin='lower'
    )

    plt.contour(
        F.T,
        levels=15,
        colors='black',
        linewidths=0.5,
        extent=extent,
        origin='lower'
    )

    plt.xlabel("tIC1",fontsize=18,fontweight='bold')
    plt.ylabel("tIC2",fontsize=18,fontweight='bold')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(-2.5,2.5)
    plt.ylim(-2.5,2.5)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)

    cbar = plt.colorbar(contour)
    cbar.set_label("Free Energy (kBT)",fontsize=15)

    plt.tight_layout()
    plt.savefig('tic1_tic2_fes_sim1_it.png',dpi=300)
    plt.show()

    # =========================
    # FEATURE CONTRIBUTION
    # =========================
    # =========================
    # Feature contribution to tIC1
    # =========================
    print("Feature contribution to tIC1 (in %)")
    weights = tica.eigenvectors[:, 0]
    contrib = np.abs(weights) / np.sum(np.abs(weights)) * 100
    idx = np.argsort(contrib)[::-1]

    for i in idx[:10]:
        print(f"{header[i+1]}: {contrib[i]:.2f}%")

    # =========================
    # Feature contribution to tIC2
    # =========================
    print("\nFeature contribution to tIC2 (in %)")
    weights = tica.eigenvectors[:, 1]
    contrib = np.abs(weights) / np.sum(np.abs(weights)) * 100
    idx = np.argsort(contrib)[::-1]

    for i in idx[:10]:
        print(f"{header[i+1]}: {contrib[i]:.2f}%")


    # =========================
    # Feature contribution to tIC3
    # =========================
    print("\nFeature contribution to tIC3 (in %)")
    weights = tica.eigenvectors[:, 2]
    contrib = np.abs(weights) / np.sum(np.abs(weights)) * 100
    idx = np.argsort(contrib)[::-1]

    for i in idx[:10]:
        print(f"{header[i+1]}: {contrib[i]:.2f}%")

    # =========================
    # Feature contribution to tIC4
    # =========================
    print("\nFeature contribution to tIC4 (in %)")
    weights = tica.eigenvectors[:, 3]
    contrib = np.abs(weights) / np.sum(np.abs(weights)) * 100
    idx = np.argsort(contrib)[::-1]

    for i in idx[:10]:
        print(f"{header[i+1]}: {contrib[i]:.2f}%")


    # =========================
    # Feature contribution to tIC5
    # =========================
    print("\nFeature contribution to tIC5 (in %)")
    weights = tica.eigenvectors[:, 4]
    contrib = np.abs(weights) / np.sum(np.abs(weights)) * 100
    idx = np.argsort(contrib)[::-1]

    for i in idx[:10]:
        print(f"{header[i+1]}: {contrib[i]:.2f}%")

if __name__ == "__main__":
    main()
