import numpy as np
import pyemma.msm as msm

# =========================
# Load tICA trajectory
# =========================

data = np.loadtxt("fes_tic1_tic2_sim1_it.dat", skiprows=1)

tic1 = data[:,0]
tic2 = data[:,1]

X = np.column_stack((tic1, tic2))

print("Total frames:", len(X))

# =========================
# Your K-means cluster centers
# =========================

centers = np.array([
[0.58099, 0.48075],
[0.86589,-0.80840],
[-1.14820,-0.06090],
[-0.34666,1.28943],
[-0.33041,0.64667],
[-2.10048,-1.73881]
])

n_clusters = len(centers)

print("Number of clusters:", n_clusters)

# =========================
# Assign frames to nearest cluster
# =========================

dtraj = np.zeros(len(X), dtype=int)

for i, point in enumerate(X):

    dist = np.sqrt(np.sum((centers - point)**2, axis=1))
    dtraj[i] = np.argmin(dist)

np.savetxt("discrete_traj.dat", dtraj, fmt="%d")

print("Discrete trajectory built.")

# =========================
# Cluster populations
# =========================

print("\nCluster populations:")

for i in range(n_clusters):

    print(f"State {i}: {np.sum(dtraj == i)} frames")

# =========================
# Build MSM (lag = 30)
# =========================

lag = 30

M = msm.estimate_markov_model([dtraj], lag=lag)

print("\nTransition matrix:")
print(M.transition_matrix)

print("\nStationary distribution:")
print(M.stationary_distribution)

print("\nActive states kept by MSM:")
print(M.active_set)

# =========================
# Map original cluster IDs → active MSM IDs
# =========================

mapping = {orig: i for i, orig in enumerate(M.active_set)}

print("\nCluster → MSM state mapping:")
print(mapping)

# =========================
# Define states (original cluster IDs)
# =========================

A_orig = [2]      # unbound
B_orig = [1]    # prebound

# convert to MSM active states
A = [mapping[x] for x in A_orig if x in mapping]
B = [mapping[x] for x in B_orig if x in mapping]

print("\nMapped A:", A)
print("Mapped B:", B)

# =========================
# Mean First Passage Times
# =========================

mfpt_AB = M.mfpt(A,B)
mfpt_BA = M.mfpt(B,A)

print("\nMFPT unbound → prebound:", mfpt_AB)
print("MFPT prebound → unbound:", mfpt_BA)

# =========================
# Transition Path Theory
# =========================

flux = msm.tpt(M, A, B)

print("\nTotal reactive flux:", flux.total_flux)

# =========================
# Dominant pathways
# =========================

paths, path_fluxes = flux.pathways(fraction=0.99)

print("\nFlux fraction      Path")
print("--------------------------------")

for i in range(len(paths)):

    frac = path_fluxes[i] / np.sum(path_fluxes)

    print(frac, paths[i])
