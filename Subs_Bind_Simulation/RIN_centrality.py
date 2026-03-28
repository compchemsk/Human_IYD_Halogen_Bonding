import MDAnalysis as mda
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
from multiprocessing import Pool
from MDAnalysis.lib.distances import distance_array

# =========================
# USER INPUT
# =========================
topology_file = "I-Tyr-binding.parm7"
trajectory_file = "sim1.nc"

DIST_CUTOFF = 6.0
ATOM_PAIR_THRESHOLD = 4
EDGE_THRESHOLD = 70.0
NPROC = 40

# =========================
# WORKER FUNCTION
# =========================
def process_frames(frame_indices):
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

    # Load once to get system info
    u = mda.Universe(topology_file, trajectory_file)
    protein = u.select_atoms("protein")
    residues = protein.residues
    n_res = len(residues)
    n_frames = len(u.trajectory)

    print(f"Total residues: {n_res}")
    print(f"Total frames: {n_frames}")
    print(f"Using {NPROC} processors")

    # Split frames
    frame_chunks = np.array_split(range(n_frames), NPROC)

    # Parallel run
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
    # PREPARE WEIGHTS (distance)
    # =========================
    for u_, v_, d in G.edges(data=True):
        d["inv_weight"] = 1.0 / d["weight"]

    # =========================
    # CENTRALITY ANALYSIS
    # =========================
    print("\nComputing centrality measures...")

    deg_cent = nx.degree_centrality(G)
    clo_cent = nx.closeness_centrality(G, distance="inv_weight")
    bet_cent = nx.betweenness_centrality(G, weight="inv_weight", normalized=True)

    print("\nResidue Centrality Measures (resid starts from 1):\n")
    print(f"{'Resid':>6} {'Resname':>8} {'Degree':>12} {'Closeness':>12} {'Betweenness':>15}")

    for node in G.nodes:
        resid = G.nodes[node]["resid"]
        resname = G.nodes[node]["resname"]

        print(f"{resid:6d} {resname:>8} {deg_cent[node]:12.6f} {clo_cent[node]:12.6f} {bet_cent[node]:15.6f}")

    # =========================
    # TOP RESIDUES
    # =========================
    def top_n(dictionary, n=10):
        return sorted(dictionary.items(), key=lambda x: x[1], reverse=True)[:n]

    print("\nTop 10 residues by Betweenness Centrality:")
    for node, val in top_n(bet_cent):
        print(f"Resid {G.nodes[node]['resid']} ({G.nodes[node]['resname']}): {val:.6f}")

    print("\nTop 10 residues by Closeness Centrality:")
    for node, val in top_n(clo_cent):
        print(f"Resid {G.nodes[node]['resid']} ({G.nodes[node]['resname']}): {val:.6f}")

    print("\nTop 10 residues by Degree Centrality:")
    for node, val in top_n(deg_cent):
        print(f"Resid {G.nodes[node]['resid']} ({G.nodes[node]['resname']}): {val:.6f}")

    # =========================
    # DOMINANT NETWORK ANALYSIS
    # =========================
    print("\nConstructing dominant interaction network (shortest paths + centrality + strong edges)...")

    dominant_nodes = []
    path = []

    if G.number_of_edges() == 0:
        print("No edges in graph!")

    else:
        # -------------------------
        # 1. Convert weight → distance
        # -------------------------
        for u_, v_, d in G.edges(data=True):
            d["dist"] = 1.0 / d["weight"]

        # -------------------------
        # 2. Shortest path (N → C)
        # -------------------------
        start_node = 0
        end_node = n_res - 1

        try:
            path = nx.shortest_path(G, source=start_node, target=end_node, weight="dist")
            print("\nShortest path (N to C):")
            print([G.nodes[n]["resid"] for n in path])
        except nx.NetworkXNoPath:
            print("\nNo path between N and C terminal!")
            path = []

        # -------------------------
        # 3. Key residues (Betweenness)
        # -------------------------
        top_k = 10
        sorted_bet = sorted(bet_cent.items(), key=lambda x: x[1], reverse=True)
        key_nodes = [node for node, _ in sorted_bet[:top_k]]

        print("\nTop Betweenness Residues:")
        print([G.nodes[n]["resid"] for n in key_nodes])

        # -------------------------
        # 4. Strong-edge subnetwork
        # -------------------------
        weights = [d["weight"] for _, _, d in G.edges(data=True)]
        threshold = np.percentile(weights, 90)  # top 10%

        G_strong = nx.Graph()
        for u, v, d in G.edges(data=True):
            if d["weight"] >= threshold:
                G_strong.add_edge(u, v, weight=d["weight"])

        print(f"\nStrong-edge threshold: {threshold:.2f}")

        # -------------------------
        # 5. Combine dominant nodes
        # -------------------------
        dominant_nodes = set()

        dominant_nodes.update(path)
        dominant_nodes.update(key_nodes)
        dominant_nodes.update(G_strong.nodes())

        dominant_nodes = list(dominant_nodes)

        print("\nCombined dominant residues:")
        print(sorted([G.nodes[n]["resid"] for n in dominant_nodes]))


    # =========================
    # VISUALIZATION
    # =========================
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)

    # Nodes
    nx.draw_networkx_nodes(G, pos, node_color="blue", node_size=200)

    # All edges (grey)
    nx.draw_networkx_edges(G, pos, edge_color="grey", width=1)

    # Highlight shortest path (black thick)
    if path and len(path) > 1:
        path_edges = list(zip(path[:-1], path[1:]))
        nx.draw_networkx_edges(
            G, pos,
            edgelist=path_edges,
            edge_color="black",
            width=3
        )

    # Highlight strong edges (slightly thicker grey/black mix)
    strong_edges = [(u, v) for u, v, d in G.edges(data=True) if d["weight"] >= threshold]
    nx.draw_networkx_edges(
        G, pos,
        edgelist=strong_edges,
        width=2
    )

    # Highlight key residues (betweenness)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=key_nodes,
        node_size=300
    )

    # Labels
    labels = {i: G.nodes[i]["resid"] for i in G.nodes}
    nx.draw_networkx_labels(G, pos, labels, font_size=8)

    plt.title("Dominant Residue Interaction Network")
    plt.axis("off")
    plt.tight_layout()
    plt.show()


    # =========================
    # HEATMAP OF RIN (80–100%)
    # =========================
    print("\nGenerating heatmap...")

    heatmap_matrix = np.zeros((n_res, n_res))

    for u, v, data in G.edges(data=True):
        weight = data["weight"]
        heatmap_matrix[u, v] = weight
        heatmap_matrix[v, u] = weight

    plt.figure(figsize=(8, 6))

    im = plt.imshow(
        heatmap_matrix,
        vmin=80, vmax=100,
        interpolation='nearest'
    )

    cbar = plt.colorbar(im)
    cbar.set_label("Contact Persistence (%)")

    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    plt.title("Residue Interaction Network Heatmap (80–100%)")

    plt.tight_layout()
    plt.show()
