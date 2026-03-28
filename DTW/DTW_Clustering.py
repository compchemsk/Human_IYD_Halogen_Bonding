import numpy as np
import matplotlib.pyplot as plt
from fastdtw import fastdtw
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.preprocessing import StandardScaler

# -------------------------------
# Load Data
# -------------------------------
rmsd_data = np.loadtxt("rmsd_full.dat")[:, 1:] #First dataset considering 100 trajectories
internal_data = np.loadtxt("int_dis_full.dat")[:, 1:] #Second dataset considering 100 trajectories

# Transpose → (num_traj, num_frames)
rmsd_data = rmsd_data.T
internal_data = internal_data.T

num_traj, num_frames = rmsd_data.shape
print(f"Loaded {num_traj} trajectories, each with {num_frames} frames.")

# -------------------------------
# Normalize Data (IMPORTANT)
# -------------------------------
scaler_rmsd = StandardScaler()
scaler_int = StandardScaler()

rmsd_data = scaler_rmsd.fit_transform(rmsd_data)
internal_data = scaler_int.fit_transform(internal_data)

# -------------------------------
# Compute DTW Distance Matrix
# -------------------------------
dtw_matrix = np.zeros((num_traj, num_traj))

for i in range(num_traj):
    for j in range(i + 1, num_traj):
        
        dist_rmsd, _ = fastdtw(rmsd_data[i], rmsd_data[j])
        dist_int, _ = fastdtw(internal_data[i], internal_data[j])
        
        # Normalize by trajectory length (important)
        dist_rmsd /= num_frames
        dist_int /= num_frames
        
        dtw_matrix[i, j] = dtw_matrix[j, i] = (dist_rmsd + dist_int) / 2

# Ensure diagonal is zero
np.fill_diagonal(dtw_matrix, 0)

# Save matrix
np.savetxt("dtw_distance_matrix.dat", dtw_matrix)

# -------------------------------
# Hierarchical Clustering
# -------------------------------
condensed = squareform(dtw_matrix)

linkage_matrix = linkage(condensed, method='ward')

# -------------------------------
# Plot Dendrogram
# -------------------------------
plt.figure(figsize=(10, 5))

dendrogram(
    linkage_matrix,
    labels=[f"Run {i+1}" for i in range(num_traj)],
    leaf_rotation=90,
    leaf_font_size=10
)

plt.title("Hierarchical Clustering of Ligand Unbinding Pathways")
plt.xlabel("Simulation Run", fontsize=18)
plt.ylabel("DTW Distance", fontsize=18)

plt.tight_layout()
plt.savefig("dtw_dendrogram.png", dpi=300)
plt.show()

# -------------------------------
# Cluster Assignment
# -------------------------------
num_clusters = 4  # adjust based on dendrogram
cluster_labels = fcluster(linkage_matrix, num_clusters, criterion='maxclust')

# Print assignments
for i, label in enumerate(cluster_labels):
    print(f"Simulation {i+1} belongs to Cluster {label}")