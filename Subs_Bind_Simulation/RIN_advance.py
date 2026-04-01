import MDAnalysis as mda
import numpy as np
import networkx as nx
from multiprocessing import Pool, cpu_count
import argparse

# -----------------------------
# STEP 1: Read files
# -----------------------------
def load_system(pdb_file, traj_file):
    u = mda.Universe(pdb_file, traj_file)
    protein = u.select_atoms("protein")

    residues = protein.residues
    n_res = len(residues)
    n_frames = len(u.trajectory)

    print(f"Number of residues: {n_res}")
    print(f"Trajectory length (frames): {n_frames}")

    return u, residues


# -----------------------------
# STEP 2: Get heavy sidechain atoms
# -----------------------------
def get_sidechain_heavy_atoms(residues):
    res_atoms = []
    for res in residues:
        atoms = res.atoms.select_atoms("not name H* and not backbone")
        res_atoms.append(atoms)
    return res_atoms


# -----------------------------
# STEP 3: Frame-wise matrix
# -----------------------------
def process_frame(args):
    frame_index, u, res_atoms = args
    u.trajectory[frame_index]

    n = len(res_atoms)
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):

            atoms_i = res_atoms[i].positions
            atoms_j = res_atoms[j].positions

            # Compute all pair distances
            diff = atoms_i[:, None, :] - atoms_j[None, :, :]
            dist = np.linalg.norm(diff, axis=2)

            # Count distances <= 6 Å
            count = np.sum(dist <= 6.0)

            if count >= 4:
                matrix[i, j] = 1
                matrix[j, i] = 1

    return matrix


# -----------------------------
# STEP 4–6: Average matrix
# -----------------------------
def compute_average_matrix(u, res_atoms, nproc):

    n_frames = len(u.trajectory)

    args = [(i, u, res_atoms) for i in range(n_frames)]

    with Pool(processes=nproc) as pool:
        matrices = pool.map(process_frame, args)

    total_matrix = np.sum(matrices, axis=0)
    avg_matrix = total_matrix / n_frames

    return avg_matrix


# -----------------------------
# STEP 7–9: Build matrices
# -----------------------------
def build_final_matrices(mat_prefinal):

    mat_final = (mat_prefinal > 0.7).astype(int)

    mat_rin = np.zeros_like(mat_prefinal)
    mat_rin[mat_final == 1] = mat_prefinal[mat_final == 1]

    return mat_final, mat_rin


# -----------------------------
# STEP 10: Build graph
# -----------------------------
def build_graph(residues, mat_final, mat_rin):

    G = nx.Graph()

    # Node positions (CA atoms)
    positions = {}

    for i, res in enumerate(residues):
        ca = res.atoms.select_atoms("name CA")[0]
        pos = ca.position
        G.add_node(i, residue=res.resid)
        positions[i] = pos[:2]  # 2D projection

    # Add edges
    n = len(residues)
    for i in range(n):
        for j in range(i + 1, n):
            if mat_final[i, j] == 1:
                weight = mat_rin[i, j]
                G.add_edge(i, j, weight=weight)

    return G, positions


# -----------------------------
# STEP 11: Centrality
# -----------------------------
def compute_centrality(G, outfile):

    bet = nx.betweenness_centrality(G, weight='weight')
    clo = nx.closeness_centrality(G)
    deg = nx.degree_centrality(G)

    with open(outfile, "w") as f:
        f.write("Residue\tBetweenness\tCloseness\tDegree\n")
        for node in G.nodes():
            f.write(f"{node}\t{bet[node]:.4f}\t{clo[node]:.4f}\t{deg[node]:.4f}\n")


# -----------------------------
# STEP 12: Save graph
# -----------------------------
def save_graph(G, outfile):
    nx.write_weighted_edgelist(G, outfile)


# -----------------------------
# MAIN
# -----------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", required=True)
    parser.add_argument("-t", "--traj", required=True)
    parser.add_argument("-n", "--nproc", type=int, default=cpu_count())

    args = parser.parse_args()

    u, residues = load_system(args.pdb, args.traj)

    res_atoms = get_sidechain_heavy_atoms(residues)

    mat_prefinal = compute_average_matrix(u, res_atoms, args.nproc)

    mat_final, mat_rin = build_final_matrices(mat_prefinal)

    G, positions = build_graph(residues, mat_final, mat_rin)

    compute_centrality(G, "centrality.txt")

    save_graph(G, "graph.edgelist")

    print("Done!")


if __name__ == "__main__":
    main()
