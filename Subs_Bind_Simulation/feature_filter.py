import numpy as np

# =========================
# INPUT / OUTPUT FILES
# =========================

input_file = "tica_msm_features.dat"
output_file = "tica_msm_features_signif.dat"

# =========================
# ALL FEATURE NAMES (order in file)
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
# SELECTED FEATURES
# =========================

selected_features = [
'RMSD_helix_A',        
'RMSD_MIA',            
'd_MIA_FRA',           
'I181_halo_binary',    
'MIA_Alpha5B',         
'd_MIA_FRB',           
'RMSD_helix_B',        
'waters_MIA_5A',       
'MIB_Alpha5A',         
'104_hydro_binary',    
'RMSD_MIB',            
'MIA_Alpha5A',         
'244*_hydro_binary',   
'248*_hydro_binary',   
'loopA_alpha5B_dist',  
'240*_hydro_binary',   
'242*_hydro_binary',   
'MIB_Alpha5B',         
'waters_helix_5A',     
'Y184_halo_binary',    
'waters_MIB_5A',       
'loopB_alpha5A_dist',  
'177_hydro_binary',    
'T171_halo_binary',    
'243*_hydro_binary',   
'waters_helix_5B',     
'258*_hydro_binary',   
'R248*_halo_binary'   
]

# =========================
# LOAD DATA
# =========================

data = np.loadtxt(input_file)

# =========================
# FIND COLUMN INDICES
# =========================

feature_indices = []

for f in selected_features:
    idx = feature_names.index(f)
    feature_indices.append(idx + 1)   # +1 because column 0 = frame

# include frame column
cols_to_keep = [0] + feature_indices

# =========================
# EXTRACT COLUMNS
# =========================

data_selected = data[:, cols_to_keep]

# =========================
# SAVE NEW FILE
# =========================

np.savetxt(output_file, data_selected, fmt="%.6f")

print("\nSaved selected features to:", output_file)
print("Columns kept:", ["frame"] + selected_features)
