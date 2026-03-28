import MDAnalysis as mda
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
from multiprocessing import Pool, cpu_count
from MDAnalysis.lib.distances import distance_array

# =========================
# USER INPUT
# =========================
topology_file = "protein.parm7"
trajectory_file = "traj.xtc"

DIST_CUTOFF = 6.0
ATOM_PAIR_THRESHOLD = 4
EDGE_THRESHOLD = 80.0
NPROC = 40  # number of processors

# =========================
# GLOBAL WORKER FUNCTION
# =========================
def process_frames(frame_indices):
    """
    Each worker loads its own Universe (IMPORTANT!)
    and processes a chunk of frames.
    """
    u = mda.Universe(topology_file, trajectory_file)
    protein = u.select_atoms("protein")
    residues = protein.residues
    n_res = len(residues)

    local_counts = np.zeros((n_res, n_res), dtype=int)

    for ts in u.trajectory[frame_indices]:

        heavy_atoms_per_res = [
            res.atoms.select_atoms("not name H*") for res in residues
        ]

        coords_per_res = [atoms.positions for atoms in heavy_atoms_per_res]

        for i, j in combinations(range(n_res), 2):

            if len(coords_per_res[i]) == 0 or len(coords_per_res[j]) == 0:
                continue

            dists = distance_array(coords_per_res[i], coords_per_res[j])
            count = np.sum(dists < DIST_CUTOFF)

            if count > ATOM_PAIR_THRESHOLD:
                local_counts[i, j] += 1
                local_counts[j, i] += 1

    return local_counts, len(frame_indices)

# =========================
# MAIN
# =========================
if __name__ == "__main__":

    u = mda.Universe(topology_file, trajectory_file)
    protein = u.select_atoms("protein")
    residues = protein.residues
    n_res = len(residues)
    n_frames = len(u.trajectory)

    print(f"Total residues: {n_res}")
    print(f"Total frames: {n_frames}")
    print(f"Using {NPROC} processors")

    # Split frames into chunks
    frame_chunks = np.array_split(range(n_frames), NPROC)

    # Parallel execution
    with Pool(processes=NPROC) as pool:
        results = pool.map(process_frames, frame_chunks)

    # Combine results
    contact_counts = np.zeros((n_res, n_res), dtype=int)
    total_frames = 0

    for counts, frames in results:
        contact_counts += counts
        total_frames += frames

    print(f"Processed frames: {total_frames}")

    # =========================
    # BUILD GRAPH
    # =========================
    G = nx.Graph()

    for i, res in enumerate(residues):
        G.add_node(i, resid=i + 1, resname=res.resname)

    for i in range(n_res):
        for j in range(i + 1, n_res):

            percent = (contact_counts[i, j] / total_frames) * 100.0

            if percent >= EDGE_THRESHOLD:
                G.add_edge(i, j, weight=percent)

    print(f"Edges after threshold: {G.number_of_edges()}")

    # =========================
    # STRONGEST PATH (N → C)
    # =========================
    start_node = 0
    end_node = n_res - 1

    for u_, v_, d in G.edges(data=True):
        d["inv_weight"] = 100.0 - d["weight"]

    try:
        path = nx.shortest_path(G, source=start_node, target=end_node, weight="inv_weight")
        print("\nDominant N → C residue path:")
        print([G.nodes[n]["resid"] for n in path])
    except nx.NetworkXNoPath:
        print("\nNo path found!")
        path = []

    # =========================
    # VISUALIZATION
    # =========================
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)

    nx.draw_networkx_nodes(G, pos, node_color="blue", node_size=200)
    nx.draw_networkx_edges(G, pos, edge_color="grey", width=1)

    if path:
        path_edges = list(zip(path[:-1], path[1:]))
        nx.draw_networkx_edges(
            G, pos,
            edgelist=path_edges,
            edge_color="black",
            width=3
        )

    labels = {i: G.nodes[i]["resid"] for i in G.nodes}
    nx.draw_networkx_labels(G, pos, labels, font_size=8)

    plt.title("Parallel Residue Interaction Network (RIN)")
    plt.axis("off")
    plt.tight_layout()
    plt.show()