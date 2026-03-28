import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.ndimage import gaussian_filter

# =========================
# Load data
# =========================

data = np.loadtxt("fes_tic1_tic2_sim1_it.dat", skiprows=1)

tic1 = data[:,0]
tic2 = data[:,1]

# =========================
# Build Free Energy Surface
# =========================

H, xedges, yedges = np.histogram2d(tic1, tic2, bins=50)

# Probability
P = H / np.sum(H)
P[P == 0] = 1e-12

# Free energy
F = -np.log(P)
F = F - np.min(F)

# Smooth the FES
F = gaussian_filter(F, sigma=1.0)

# =========================
# Assign FES value to each trajectory point
# =========================

x_bin = np.digitize(tic1, xedges) - 1
y_bin = np.digitize(tic2, yedges) - 1

# Keep indices inside range
x_bin = np.clip(x_bin, 0, F.shape[0]-1)
y_bin = np.clip(y_bin, 0, F.shape[1]-1)

# Free energy of each point
F_point = F[x_bin, y_bin]

# =========================
# Apply free energy cutoff
# =========================

mask = F_point < 10

tic1_cut = tic1[mask]
tic2_cut = tic2[mask]

X = np.column_stack((tic1_cut, tic2_cut))

# =========================
# K-means clustering
# =========================

n_clusters = 6

kmeans = KMeans(n_clusters=n_clusters, random_state=0)

labels = kmeans.fit_predict(X)

centers = kmeans.cluster_centers_

# =========================
# Save cluster centers
# =========================

with open("cluster_centers.dat","w") as f:

    f.write("Cluster   tIC1   tIC2\n")

    for i,(cx,cy) in enumerate(centers):

        line = f"{i}   {cx:.5f}   {cy:.5f}"

        print(line)
        f.write(line+"\n")

# =========================
# Plot Free Energy Surface
# =========================

plt.figure(figsize=(8,6))

sc = plt.scatter(
    tic1,
    tic2,
    c=F_point,
    cmap="jet",
    s=20
)

# Plot cluster centers
plt.scatter(
    centers[:,0],
    centers[:,1],
    color="white",
    s=200,
    edgecolors="black",
    marker="o",
    label="Cluster centers"
)

plt.xlabel("tIC1", fontsize=18, fontweight='bold')
plt.ylabel("tIC2", fontsize=18, fontweight='bold')

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

ax = plt.gca()

for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(2)

cbar = plt.colorbar(sc)
cbar.set_label("Free Energy (kBT)", fontsize=15)

plt.legend()

plt.tight_layout()
plt.savefig('tic1_tic2_kmeans_sim1_it.png',dpi=300)
plt.show()

# =========================
# Find 100 nearest points to each center
# =========================

n_nearest = 100

for i, center in enumerate(centers):

    # compute distance from all points
    dist = np.sqrt((tic1 - center[0])**2 + (tic2 - center[1])**2)

    # get indices of 100 closest points
    nearest_idx = np.argsort(dist)[:n_nearest]

    # save the points
    fname = f"cluster_{i}_nearest100.dat"

    with open(fname, "w") as f:
        f.write("Index   tIC1   tIC2   Distance\n")

        for idx in nearest_idx:
            f.write(f"{idx}   {tic1[idx]:.5f}   {tic2[idx]:.5f}   {dist[idx]:.5f}\n")

    print(f"Saved {fname}")
