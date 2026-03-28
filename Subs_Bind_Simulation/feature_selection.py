import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array
from multiprocessing import Pool
from tqdm import tqdm

# =========================
# INPUT FILES
# =========================

topology = "I-Tyr-binding.parm7"
trajectories = ["sim1.nc"]

# =========================
# PRECOMPUTE HBONDS (SERIAL)
# =========================

u = mda.Universe(topology, trajectories)

align.AlignTraj(
    u,
    u,
    select="protein and backbone",
    in_memory=True
).run()

hbond_analysis = HydrogenBondAnalysis(
    universe=u,
    donors_sel="protein or resname MIA or resname MIB",
    hydrogens_sel="name H*",
    acceptors_sel="protein or resname MIA or resname MIB",
    d_h_cutoff=1.2,
    d_a_cutoff=3.0,
    d_h_a_angle_cutoff=150,
    update_selections=False
)

print("Starting hydrogen bond analysis...")
hbond_analysis.run()

hbonds_by_frame = {}

for hbond in hbond_analysis.results.hbonds:
    frame = int(hbond[0])
    hbonds_by_frame.setdefault(frame, []).append(hbond)
print("Hydrogen bond analysis finished.")
print("Starting parallel feature extraction...")


def process_frame(frame):

    u = mda.Universe(topology, trajectories)
    u.trajectory[frame]

    protein = u.select_atoms("protein")
    waters  = u.select_atoms("resname WAT")

    MIA = u.select_atoms("resname MIA")
    MIB = u.select_atoms("resname MIB")

    MIA_ring = u.select_atoms("resname MIA and name C4 C7 C8 C5 C6 C9")
    MIB_ring = u.select_atoms("resname MIB and name C4 C7 C8 C5 C6 C9")

    FRA_core = u.select_atoms("resname FRA and name N1 C2 N3 C4 C4a C10a")
    FRB_core = u.select_atoms("resname FRB and name N1 C2 N3 C4 C4a C10a")

    helix_A = u.select_atoms("resid 73-118 and backbone")
    helix_B = u.select_atoms("resid 223-268 and backbone")

    loop_B   = u.select_atoms("resid 278-292")
    alpha5_A = u.select_atoms("resid 97-107")

    loop_A   = u.select_atoms("resid 128-142")
    alpha5_B = u.select_atoms("resid 247-257")

    data = [frame]

    # ----------------------
    # COM distances
    # ----------------------
    MIA_com = MIA_ring.center_of_mass()
    MIB_com = MIB_ring.center_of_mass()
    FRA_com = FRA_core.center_of_mass()
    FRB_com = FRB_core.center_of_mass()

    data += [
        np.linalg.norm(MIA_com - FRA_com),
        np.linalg.norm(MIA_com - FRB_com),
        np.linalg.norm(MIB_com - FRA_com),
        np.linalg.norm(MIB_com - FRB_com)
    ]

    # ----------------------
    # Minimum distances
    # ----------------------
    data += [
        distance_array(MIA.positions, alpha5_A.positions).min(),
        distance_array(MIB.positions, alpha5_B.positions).min(),
        distance_array(MIA.positions, alpha5_B.positions).min(),
        distance_array(MIB.positions, alpha5_A.positions).min()
    ]

    # ----------------------
    # Water counting
    # ----------------------
    cutoff = 3.5
    def count_waters(site_atoms):
        dist = distance_array(site_atoms.positions, waters.positions)
        close_atoms = waters[np.any(dist < cutoff, axis=0)]
        return len(np.unique(close_atoms.resids))

    data += [
        count_waters(MIA),
        count_waters(MIB),
        count_waters(alpha5_A),
        count_waters(alpha5_B)
    ]

    # ----------------------
    # RMSD
    # ----------------------
    u.trajectory[0]
    ref_A = alpha5_A.positions.copy()
    ref_B = alpha5_B.positions.copy()
    ref_A1 = MIA.positions.copy()
    ref_B1 = MIB.positions.copy()
    u.trajectory[frame]

    data += [
        rmsd(alpha5_A.positions, ref_A, superposition=False),
        rmsd(alpha5_B.positions, ref_B, superposition=False),
        rmsd(MIA.positions, ref_A1, superposition=False),
        rmsd(MIB.positions, ref_B1, superposition=False),
    ]

    # ----------------------
    # Loop distances
    # ----------------------
    data += [
        np.linalg.norm(loop_B.center_of_mass() - alpha5_A.center_of_mass()),
        np.linalg.norm(loop_A.center_of_mass() - alpha5_B.center_of_mass())
    ]

    # ----------------------
    # Hydrogen bonds
    # ----------------------
    def detect_hbond(frame,resid,lig):
        if frame not in hbonds_by_frame:
            return 0

        residue_atoms = u.select_atoms(f"protein and resid {resid}")
        ligand_atoms = u.select_atoms(f"resname {lig}")

        for h in hbonds_by_frame[frame]:
            donor = int(h[1])
            acceptor = int(h[3])
            if ((donor in residue_atoms.indices and acceptor in ligand_atoms.indices) or
                (donor in ligand_atoms.indices and acceptor in residue_atoms.indices)):
                return 1
        return 0

    data += [
        detect_hbond(frame,390,"MIA"),
        detect_hbond(frame,392,"MIA"),
        detect_hbond(frame,393,"MIA"),
        detect_hbond(frame,394,"MIA"),
        detect_hbond(frame,398,"MIA"),
        detect_hbond(frame,408,"MIA"),
        detect_hbond(frame,107,"MIB"),
        detect_hbond(frame,34,"MIB"),
        detect_hbond(frame,356,"MIB"),
        detect_hbond(frame,358,"MIB"),
        detect_hbond(frame,98,"MIB")
    ]

    return data


# =========================
# PARALLEL EXECUTION
# =========================

if __name__ == "__main__":

    n_frames = len(u.trajectory)

    with Pool(40) as pool:

        results = list(
            tqdm(
                pool.imap(process_frame, range(n_frames)),
                total=n_frames
            )
        )

    results.sort(key=lambda x: x[0])

    outfile = open("tica_msm_features_new.dat","w")

    for row in results:
        outfile.write(" ".join(map(str,row))+"\n")

    outfile.close()


import mdtraj as md
import numpy as np

def calculate_angle(vec1, vec2):
    """Calculate the angle (in degrees) between two vectors."""
    cos_theta = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

# Load trajectory and topology
traj = md.load('sim1.nc', top='4ttc-unbiased-fmnh2-ityr.parm7')
top = traj.topology

# Distance cutoffs in nm (mdtraj positions are in nm)
distance_criteria = {'N': 0.397, 'O': 0.382, 'S': 0.415, 'aromatic': 0.400}
angle_range = (140, 180)  # degrees

# Residues of interest
residues_of_interest = {
    "MIA": [397, 107, 327, 408],
    "MIB": [110, 100, 204, 113]
}

# Ligand halogen and carbon indices
ligand_indices = {
    "MIA": {"I": 7376, "C": 7361},
    "MIB": {"I": 7400, "C": 7385}
}

# Define aromatic rings
ring_atoms = {
    'TYR': ['CD1','CD2','CE1','CE2','CG','CZ'],
    'PHE': ['CD1','CD2','CE1','CE2','CG','CZ'],
    'TRP': ['CZ2','CH2','CZ3','CE3','CD2','CE2','CG','CD1','NE1']
}

# Precompute donor atom indices per residue
halogen_donors = {}
for residue in top.residues:
    donor_list = []
    # Aromatic rings
    if residue.name in ring_atoms:
        aromatic_indices = [atom.index for atom in residue.atoms if atom.name in ring_atoms[residue.name]]
        if aromatic_indices:
            donor_list.append(("aromatic", aromatic_indices))
    # N, O, S atoms
    for atom in residue.atoms:
        if atom.element.symbol == 'N':
            donor_list.append(('N', [atom.index]))
        elif atom.element.symbol == 'O':
            donor_list.append(('O', [atom.index]))
        elif atom.element.symbol == 'S':
            donor_list.append(('S', [atom.index]))
    if donor_list:
        halogen_donors[residue.index] = donor_list

# Initialize results dictionary
data = {resid: [] for reslist in residues_of_interest.values() for resid in reslist}

# Iterate through all frames
for frame_idx in range(traj.n_frames):
    xyz = traj.xyz[frame_idx]
    
    for ligand_type, residues in residues_of_interest.items():
        iodine_pos = xyz[ligand_indices[ligand_type]["I"]]
        carbon_pos = xyz[ligand_indices[ligand_type]["C"]]
        
        for resid in residues:
            hb_found = 0
            if resid in halogen_donors:
                for donor_type, donor_indices in halogen_donors[resid]:
                    if donor_type == 'aromatic':
                        donor_pos = np.mean(xyz[donor_indices], axis=0)
                    else:
                        donor_pos = xyz[donor_indices[0]]
                    
                    distance = np.linalg.norm(iodine_pos - donor_pos)
                    if distance <= distance_criteria[donor_type]:
                        vec1 = carbon_pos - iodine_pos
                        vec2 = donor_pos - iodine_pos
                        angle = calculate_angle(vec1, vec2)
                        if angle_range[0] <= angle <= angle_range[1]:
                            hb_found = 1
                            break  # Stop after first detected halogen bond
            data[resid].append(hb_found)

# Order results same as your original call
result = [
    data[397], data[107], data[327], data[408],  # MIA
    data[110], data[100], data[204], data[113]   # MIB
]

# Convert result to 2D array (rows=frames, columns=residues)
framewise_array = np.array(result).T  # transpose

# File path to save
output_file = "halogen_bonds.txt"

# Save to file without header
np.savetxt(output_file, framewise_array, fmt='%d', delimiter='\t')


import numpy as np
import os

# Load the two files
tica_data = np.loadtxt("tica_msm_features_new.dat")  # shape: (frames, features)
halogen_data = np.loadtxt("halogen_bonds.txt")       # shape: (frames, residues)

# Concatenate side by side (columns)
combined_data = np.hstack((tica_data, halogen_data))

# Save to new file
np.savetxt("tica_msm_features.dat", combined_data, fmt='%.6f', delimiter='\t')

# Remove the old files
os.remove("tica_msm_features_new.dat")
os.remove("halogen_bonds.txt")

print("Feature Extraction Complete.")
