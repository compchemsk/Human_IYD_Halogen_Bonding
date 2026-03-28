import numpy as np
import pyemma
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

# =========================
# PARAMETERS
# =========================

feature_file = "tica_msm_features.dat"

lag = 30
dim = 5
chunk_size = 1000
equilibration = 2420
cv_splits = 10

np.random.seed(42)

# =========================
# FEATURE NAMES
# =========================

feature_names = [
"d_MIA_FRA","d_MIA_FRB",
"d_MIB_FRA","d_MIB_FRB",
"MIA_Alpha5A","MIB_Alpha5B","MIA_Alpha5B","MIB_Alpha5A",
"waters_MIA_5A","waters_MIB_5A","waters_helix_5A","waters_helix_5B",
"RMSD_helix_A","RMSD_helix_B","RMSD_MIA","RMSD_MIB",
"loopB_alpha5A_dist","loopA_alpha5B_dist",
"240*_hydro_binary","242*_hydro_binary","243*_hydro_binary",
"244*_hydro_binary","248*_hydro_binary","258*_hydro_binary",
"177_hydro_binary","104_hydro_binary","206*_hydro_binary",
"208*_hydro_binary","168_hydro_binary",
"R248*_halo_binary","T178_halo_binary","T178*_halo_binary",
"L259*_halo_binary","I181_halo_binary","T171_halo_binary",
"P275_halo_binary","Y184_halo_binary"
]

# =========================
# LOAD DATA
# =========================

print("Loading data...")

data = np.loadtxt(feature_file)
data = data[equilibration:]

X = data[:,1:38]

print("Feature matrix:", X.shape)

# =========================
# CREATE PSEUDO TRAJECTORIES
# =========================

trajectories = []

for i in range(0, len(X), chunk_size):

    chunk = X[i:i+chunk_size]

    if len(chunk) > 2*lag:
        trajectories.append(chunk)

print("Pseudo trajectories:", len(trajectories))

# =========================
# SCALE FEATURES
# =========================

scaled_trajs = []

for traj in trajectories:

    scaler = StandardScaler()

    scaled_trajs.append(scaler.fit_transform(traj))

# =========================
# CROSS VALIDATED VAMP2
# =========================

def score_cv(data, dim, lag):

    scores = []

    for _ in range(cv_splits):

        idx = np.random.permutation(len(data))
        split = len(data)//2

        train = [data[i] for i in idx[:split]]
        test  = [data[i] for i in idx[split:]]

        vamp = pyemma.coordinates.vamp(train, lag=lag, dim=dim)

        scores.append(vamp.score(test))

    return np.mean(scores), np.std(scores)

# =========================
# INDIVIDUAL FEATURE VAMP2
# =========================

print("\nComputing individual feature scores\n")

feature_scores = []
feature_errors = []

for i,fname in enumerate(feature_names):

    feature_trajs = []

    for traj in scaled_trajs:

        f = traj[:,i].reshape(-1,1)

        if np.std(f) < 1e-8:
            continue

        feature_trajs.append(f)

    if len(feature_trajs) < 2:

        feature_scores.append(0)
        feature_errors.append(0)

        print(fname,"skipped")
        continue

    mean,std = score_cv(feature_trajs, dim=1, lag=lag)

    feature_scores.append(mean)
    feature_errors.append(std)

    print(f"{fname:25s} mean={mean:.4f}")

# =========================
# BAR PLOT (INDIVIDUAL FEATURES)
# =========================

order = np.argsort(feature_scores)[::-1]

sorted_names = [feature_names[i] for i in order]
sorted_scores = [feature_scores[i] for i in order]
sorted_errors = [feature_errors[i] for i in order]

plt.figure(figsize=(14,6))

plt.bar(sorted_names, sorted_scores, yerr=sorted_errors)

plt.xticks(rotation=90,fontsize=15,fontweight='bold')
plt.yticks(fontsize=15,fontweight='bold')

plt.ylabel("VAMP2 Score",fontsize=20,fontweight='bold')
# Axis styling
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(2)


plt.tight_layout()

plt.savefig("sim1_vamp2_indiv.png", dpi=300)

plt.show()

# =========================
# FEATURE SELECTION USING CORRELATION FILTERING
# =========================

top_n = 30
corr_threshold = 0.8

print("\nSelecting top features based on VAMP2 score\n")

# indices of top N VAMP2 features
top_indices = order[:top_n]

top_features = [feature_names[i] for i in top_indices]

print("Top features from VAMP2 ranking:\n")

for i,f in enumerate(top_features):
    print(i+1,f)

# =========================
# BUILD MATRIX OF THESE FEATURES
# =========================

X_top = X[:, top_indices]

# correlation matrix
#corr_matrix = np.corrcoef(X_top.T)

from sklearn.feature_selection import mutual_info_regression

n_feat = X_top.shape[1]

# =========================
# MUTUAL INFORMATION MATRIX
# =========================

corr_matrix = np.zeros((n_feat, n_feat))

for i in range(n_feat):

    y = X_top[:, i]

    corr_matrix[:, i] = mutual_info_regression(
        X_top,
        y,
        discrete_features=False
    )

# make symmetric
coor_matrix = (corr_matrix + corr_matrix.T) / 2

# =========================
# REMOVE HIGHLY CORRELATED FEATURES
# =========================

selected_indices = []
removed_indices = []

for i in range(len(top_indices)):

    keep = True

    for j in selected_indices:

        if abs(corr_matrix[i,j]) > corr_threshold:
            keep = False
            removed_indices.append(i)
            break

    if keep:
        selected_indices.append(i)

final_indices = [top_indices[i] for i in selected_indices]

final_features = [feature_names[i] for i in final_indices]

print("\nFinal selected feature set (after correlation filtering):\n")

for i,f in enumerate(final_features):
    print(i+1,f)

# =========================
# CORRELATION HEATMAP
# =========================

import matplotlib.pyplot as plt

plt.figure(figsize=(10,8))

plt.imshow(corr_matrix, cmap='seismic', vmin=0, vmax=1)

plt.colorbar(label="Correlation")

plt.xticks(range(len(top_features)), top_features, rotation=90, fontsize=12,fontweight='bold')
plt.yticks(range(len(top_features)), top_features, fontsize=12,fontweight='bold')

plt.tight_layout()

plt.savefig("sim1_feat_coor_vamp2_it.png", dpi=300)

plt.show()

# =========================
# FINAL FEATURE IMPORTANCE PLOT
# =========================

final_scores = [feature_scores[i] for i in final_indices]
final_errors = [feature_errors[i] for i in final_indices]

plt.figure(figsize=(10,6))

plt.bar(final_features, final_scores, yerr=final_errors)

plt.xticks(rotation=90, fontsize=15, fontweight='bold')
plt.yticks(fontsize=15, fontweight='bold')

plt.ylabel("VAMP2 Score", fontsize=20, fontweight='bold')

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(2)

plt.tight_layout()

plt.savefig("sim1_feat_it.png", dpi=300)

plt.show()

# =========================
# PRINT FINAL ORDERING
# =========================

print("\nFinal feature ordering (best → worst among selected):\n")

sorted_final = sorted(
    zip(final_features, final_scores),
    key=lambda x: x[1],
    reverse=True
)

for i,(f,s) in enumerate(sorted_final):
    print(f"{i+1:2d}  {f:25s}  VAMP2={s:.4f}")

