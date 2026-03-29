# 🧬 Human_IYD_Halogen_Bonding

Human IYD shows **no halogen bonding in crystal structures**, yet kinetics reveal **strong halogen-dependent binding**. Heavier halotyrosines bind tighter and associate faster than lighter halotyrosines, while *k_off* changes little—suggesting halogens mainly boost *k_on* via **transient interactions**. Static structures miss these dynamic effects.  
📌 Citation: S. Karmakar and S. Mishra, J. Chem. Inf. Model., 2026 (under review) (Cite: https://doi.org/10.5281/zenodo.19310318)

---

## 🔬 Study Overview

Previous computational and experimental studies showed **absence of halogen bonding** in the active site. However, kinetic trends suggest a **hidden role during substrate recognition**.

### 🧪 Systems Studied:
- **I-Tyr** → Strong halogen + hydrogen bonding  
- **Cl-Tyr** → Weak halogen + strong hydrogen bonding  
- **I-Phenol** → Strong halogen + weak hydrogen bonding  

---

# ⚙️ A. Parameter Preparation
- Follow protocol from:  
  https://github.com/compchemsk/Human_IYD_Dynamics.git  
- Parameter files are available in **Parameters/**
- MD input files same as previous repository

---

# 🧭 B. Classical MD & iSMD Setup

## 🌀 Step 1: Classical MD
- Run **1 μs MD simulation** for each system  
- Extract **minimum free energy structure** using PCA (PC1 vs PC2)

## 🧩 Step 2: Pre-equilibration for iSMD
- Input files → iSMD/
- Run standard MD commands

Output:
Pre-Equilibration.rst7

## 🚀 Step 3: Pulling Simulations (iSMD)

Force constants:
0.1 → 4 kcal·mol⁻¹·Å⁻²

Run:
chmod +x Pull_jobs.sh  
./Pull_jobs.sh

---

# 🎯 C. Goldilocks Force Constant Selection

Directory structure:
k0.1 → k4

## 📊 Work Distribution
chmod +x Distribution.sh  
./Distribution.sh  
python work-dist-plot.py  

## 📉 Free Energy Distribution
python free_energy_smd.py  

## 💪 Force Distribution
(loop over k values and run k_F.sh, then python k_F_plot.py)

## 🧮 PMF Distribution
(loop over k values and run k_PMF.sh, then python k_PMF_plot.py)

## 🔗 Internal Coordinates
(loop and run int-coord.sh, then int-coord-fes.sh)

## 🔄 Sequence of Events
chmod +x seq-of-events.sh  
./seq-of-events.sh  

---

# 📉 D. Kinetics via DHS Model
Inputs:
- Rupture force
- Force SD
- Unbinding time
- Time SD

---

# 🧭 E. Pathway Identification (DTW)

chmod +x int-dist-full.sh  
./int-dist-full.sh  

chmod +x rmsd_full.sh  
./rmsd_full.sh  

pip install fastdtw scipy seaborn matplotlib numpy  
python DTW_Clustering.py  

---

# 🪟 F. Umbrella Sampling

Optimal:
- Early: k=10, spacing=0.2 Å  
- Late: k=2.5, spacing=1 Å  

25 ns per window

Run:
target.sh  
rst7-input.sh  
name-change.sh  
job-us-run.sh  

PMF:
./npwham-run.sh  
python pmf-plot.py  

---

# 🧪 G. Halogen Bond Analysis

Criteria:
- Distance cutoff  
- Angle 140°–180°  

Run:
./halogen-bond.sh  

Additional:
halo-amino.sh  
halo-count.sh  
halo-plot-data.py  
halo-along-pathway.py  

QM Validation includes NBO, DFT, ALMO-EDA.

---

# 💧 H. Hydrogen Bond Analysis

Criteria:
- Donor–H ≤ 1.2 Å  
- Donor–Acceptor ≤ 3.0 Å  
- Angle ≥ 150°  

---

# 🔥 I. Free Energy Contributions
- Solvation  
- Halogen bonding  
- Hydrogen bonding  
- Sterics  

---

# 🧬 J. Substrate Binding Simulations

- Substrates placed 50 Å away  

## 🌐 RIN
- Nodes: residues  
- Edges: heavy atom contacts  
- Persistence ≥ 80%  

## 🧠 TICA Workflow
feature_selection.py → lagfind.py → VAMP.py → feature_filter.py → lagfind_filter.py → MSM.py → kmeans.py → CK.py → TPT.py  

---

# Final Note
Please cite appropriately if you use this workflow.
