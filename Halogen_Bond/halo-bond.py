import mdtraj as md
import numpy as np

def calculate_angle(vec1, vec2):
    """Calculate the angle (in degrees) between two vectors."""
    cos_theta = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

# Load trajectory and topology
traj = md.load('xxxx-us-10ns.nc', top='hIYD-I-Tyr.parm7')
top = traj.topology

# Define indices
iodine_index = 7324
carbon_index = 7309

# Convert distance criteria from Å to nm
distance_criteria = {'N': 0.397, 'O': 0.382, 'S': 0.415, 'aromatic': 0.400} #These values are specific for iodine, will change in case of chlorine
angle_range = (140, 180)

# Identify potential halogen bond donors
halogen_donors = []
ring_center_atoms = {'TYR': [['CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ']],
                     'PHE': [['CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ']],
                     'TRP': [['CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CE2'],
                             ['CE2', 'CD2', 'CG', 'CD1', 'NE1']]}  # Two ring centers for TRP

# Identify potential donors in the system
for residue in top.residues:
    if residue.name in ring_center_atoms:
        for atom_group in ring_center_atoms[residue.name]:
            atom_indices = [atom.index for atom in residue.atoms if atom.name in atom_group]
            if atom_indices:
                halogen_donors.append((atom_indices, 'aromatic'))
    for atom in residue.atoms:
        if atom.name.startswith('N'):
            halogen_donors.append((atom.index, 'N'))
        elif atom.name.startswith('O'):
            halogen_donors.append((atom.index, 'O'))
        elif atom.name.startswith('S'):
            halogen_donors.append((atom.index, 'S'))

# Store results
halogen_bonds = {}

# Iterate through trajectory frames (0 to 350, i.e., 351 frames)
for frame_idx in range(500):
    xyz = traj.xyz[frame_idx]
    iodine_pos = xyz[iodine_index]
    carbon_pos = xyz[carbon_index]

    for donor_info, donor_type in halogen_donors:
        if donor_type == 'aromatic':
            donor_positions = xyz[donor_info]
            donor_pos = np.mean(donor_positions, axis=0)  # Center of mass of the ring
            atom_type_info = "aromatic ring"  # Special case for aromatic systems
        else:
            donor_pos = xyz[donor_info]
            atom_type_info = top.atom(donor_info).name  # Get atom name (e.g., ND1, OE1, etc.)

        # Calculate distance between iodine and donor
        distance = np.linalg.norm(iodine_pos - donor_pos)

        # Check if distance is within the defined criteria
        if distance <= distance_criteria[donor_type]:
            vec1 = carbon_pos - iodine_pos  # Carbon to Iodine vector
            vec2 = donor_pos - iodine_pos  # Donor to Iodine vector
            angle = calculate_angle(vec1, vec2)

            # Check if angle falls within the acceptable range
            if angle_range[0] <= angle <= angle_range[1]:
                if donor_type == 'aromatic':
                    residue_info = f"{top.atom(donor_info[0]).residue} ({donor_type}, {atom_type_info})"
                else:
                    residue_info = f"{top.atom(donor_info).residue} ({donor_type}, {atom_type_info})"
                if frame_idx not in halogen_bonds:
                    halogen_bonds[frame_idx] = []
                halogen_bonds[frame_idx].append((residue_info, distance, donor_pos, angle))

# Print results
if halogen_bonds:
    print("For trajectory xxxx")
    for frame, residues in sorted(halogen_bonds.items()):
        print(f"Frame {frame}:")
        for res_info, dist, d_pos, ang in residues:
            print(f"  {res_info}, Distance: {dist:.3f} nm, Donor Pos: {d_pos}, Angle: {ang:.2f}°")
else:
    print("For trajectory xxxx")
    print("No halogen bonds detected in this trajectory")
