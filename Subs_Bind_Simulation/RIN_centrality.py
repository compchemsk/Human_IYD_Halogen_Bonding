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
topology_file = "I-Tyr-binding.parm7" #provide your parameter file
trajectory_file = "sim1.nc" #provide your trajectory file

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
    print("\nConstructing dominant interaction network (robust pathway)...")

    dominant_nodes = set()
    dominant_edges = set()

    if G.number_of_edges() == 0:
        print("No edges in graph!")

    else:
        # -------------------------
        # 1. Convert weight → distance
        # -------------------------
        for u, v, d in G.edges(data=True):
            d["dist"] = 1.0 / (d["weight"] + 1e-6)  # avoid divide-by-zero

        # -------------------------
        # 2. k-shortest paths (robust)
        # -------------------------
        start_node = 0
        end_node = n_res - 1
        k_paths = 5
        paths = []

        try:
            paths_gen = nx.shortest_simple_paths(G, start_node, end_node, weight="dist")
            for i, p in enumerate(paths_gen):
                if i >= k_paths:
                    break
                paths.append(p)

            print(f"\nTop {len(paths)} shortest paths:")
            for p in paths:
                print([G.nodes[n]["resid"] for n in p])

        except nx.NetworkXNoPath:
            print("\nNo path between N and C terminal!")
            paths = []

        # Collect path nodes & edges
        for p in paths:
            dominant_nodes.update(p)
            dominant_edges.update(zip(p[:-1], p[1:]))

        # -------------------------
        # 3. Node Betweenness
        # -------------------------
        bet_cent = nx.betweenness_centrality(G, weight="dist", normalized=True)

        top_k = max(5, int(0.1 * G.number_of_nodes()))
        sorted_bet = sorted(bet_cent.items(), key=lambda x: x[1], reverse=True)
        key_nodes = [node for node, _ in sorted_bet[:top_k]]

        print("\nTop Betweenness Residues:")
        print([G.nodes[n]["resid"] for n in key_nodes])

        dominant_nodes.update(key_nodes)

        # -------------------------
        # 4. Edge Betweenness (VERY IMPORTANT)
        # -------------------------
        edge_bet = nx.edge_betweenness_centrality(G, weight="dist")

        edge_thresh = np.percentile(list(edge_bet.values()), 90)

        key_edges = [e for e, val in edge_bet.items() if val >= edge_thresh]

        print(f"\nEdge betweenness threshold: {edge_thresh:.4f}")

        dominant_edges.update(key_edges)

        # -------------------------
        # 5. Strong-edge subnetwork
        # -------------------------
        weights = [d["weight"] for _, _, d in G.edges(data=True)]
        weight_thresh = np.percentile(weights, 90)

        strong_edges = [(u, v) for u, v, d in G.edges(data=True)
                        if d["weight"] >= weight_thresh]

        print(f"Strong-edge threshold: {weight_thresh:.2f}")

    # =========================
    # VISUALIZATION
    # =========================
    plt.figure(figsize=(12, 9))

    pos = nx.spring_layout(G, seed=42)

    # --- Base network ---
    nx.draw_networkx_edges(G, pos, edge_color="lightgrey", width=0.8)
    nx.draw_networkx_nodes(G, pos, node_color="lightblue", node_size=150)

    # --- Dominant edges (red) ---
    nx.draw_networkx_edges(
        G, pos,
        edgelist=list(dominant_edges),
        width=3
    )

    # --- Dominant nodes (larger) ---
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=list(dominant_nodes),
        node_size=400
    )

    # --- Labels ---
    labels = {i: G.nodes[i]["resid"] for i in dominant_nodes}
    nx.draw_networkx_labels(G, pos, labels, font_size=8)

    plt.title("Dominant Interaction Pathway (Robust Network)")
    plt.axis("off")
    plt.tight_layout()
    plt.show()
