<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">

</head>

<body>

<div class="container">

<h1>Human_IYD_Halogen_Bonding</h1>


<p>
Human IYD shows no halogen bonding in crystal structures, yet kinetics reveal strong halogen-dependent binding. 
Heavier halotyrosines bind tighter and associate faster than lighter halotyrosines, while <i>koff</i> changes little, 
suggesting halogens mainly boost <i>kon</i> via transient interactions.
</p>

<p class="highlight">
Citation: S. Karmakar and S. Mishra, J. Chem. Inf. Model., 2026 (under review)
</p>

<h2>Background</h2>
<p>
Previous studies showed no halogen bonding in static structures. However, kinetics indicate a hidden role during substrate recognition.
</p>

<ul>
<li>I-Tyr → Strong halogen + strong hydrogen bonding</li>
<li>Cl-Tyr → Weak halogen + strong hydrogen bonding</li>
<li>I-Phenol → Strong halogen + weak hydrogen bonding</li>
</ul>

<h2>A. System Preparation</h2>
<p>
Parameter files were created using the Human_IYD_Dynamics repository.
</p>

<h2>B. Classical MD and iSMD Setup</h2>

<ul>
<li>Run 1 μs MD simulation</li>
<li>Extract minimum free energy structure</li>
<li>Perform pre-equilibration</li>
</ul>

<p><b>Output:</b></p>
<pre>Pre-Equilibration.rst7</pre>

<p>This file is used for pulling simulations.</p>

<pre>
chmod +x Pull_jobs.sh
./Pull_jobs.sh
</pre>

<h2>C. Goldilocks Force Constant Analysis</h2>

<h3>1. Work Distribution</h3>
<pre>
chmod +x Distribution.sh
./Distribution.sh
python work-dist-plot.py
</pre>

<h3>2. Free Energy Distribution</h3>
<pre>python free_energy_smd.py</pre>

<h3>3. Force Distribution</h3>
<pre>
for np in {0.1,0.25,0.5,1,2,3,4}
...
python k_F_plot.py
</pre>

<h3>4. PMF Distribution</h3>
<pre>
for np in {0.1,0.25,0.5,1,2,3,4}
...
python k_PMF_plot.py
</pre>

<h3>5. Internal Degrees of Freedom</h3>
<pre>
for np in {0.1,0.25,0.5,1,2,3,4}
...
</pre>

<h3>6. Sequence of Events</h3>
<pre>
chmod +x seq-of-events.sh
./seq-of-events.sh
</pre>

<h2>D. DHS Kinetics</h2>
<ul>
<li>Rupture forces</li>
<li>Standard deviations</li>
<li>Unbinding time</li>
</ul>

<h2>E. DTW Clustering</h2>

<pre>
chmod +x int-dist-full.sh
./int-dist-full.sh

chmod +x rmsd_full.sh
./rmsd_full.sh

pip install fastdtw scipy seaborn matplotlib numpy
python DTW_Clustering.py
</pre>

<p>
Clustering based on RMSD and COM distance using Ward linkage.
</p>

<h2>F. Umbrella Sampling</h2>

<ul>
<li>k = 10 kcal/mol (0.2 Å spacing)</li>
<li>k = 2.5 kcal/mol (1 Å spacing)</li>
<li>25 ns/window</li>
</ul>

<pre>
cpptraj → dis.dat
./target.sh
./rst7-input.sh
./name-change.sh
./job-us-run.sh
</pre>

<pre>
./npwham-run.sh
python pmf-plot.py
</pre>

<h2>G. Halogen Bond Analysis</h2>

<p>Criteria: distance + angle (140°–180°)</p>

<pre>
./halogen-bond.sh → halo.out
./halo-amino.sh
./halo-count.sh
python halo-plot-data.py
python halo-along-pathway.py
</pre>

<h3>QM Validation</h3>
<ul>
<li>DFT + NBO</li>
<li>E(2) ≥ 0.5 kcal/mol</li>
</ul>

<h2>H. Hydrogen Bond Analysis</h2>
<ul>
<li>D–H ≤ 1.2 Å</li>
<li>D–A ≤ 3.0 Å</li>
<li>Angle ≥ 150°</li>
</ul>

<h2>I. Free Energy Contributions</h2>
<ul>
<li>Solvation</li>
<li>Halogen bonding</li>
<li>Hydrogen bonding</li>
<li>Steric crowding</li>
</ul>

<h2>J. Substrate Binding Simulation</h2>

<p>Substrates placed 50 Å away.</p>

<h3>Residue Interaction Network</h3>
<ul>
<li>>4 heavy atom contacts within 6 Å</li>
<li>>80% persistence</li>
</ul>

<h3>TICA + MSM Workflow</h3>
<pre>
feature_selection.py → lagfind.py → VAMP.py → feature_filter.py → lagfind_filter.py → MSM.py → kmeans.py → CK.py → TPT.py
</pre>

<p>
Transition network identifies pathway from unbound to pre-bound state.
</p>

<div class="footer">
Thank you. Please cite appropriately.
</div>

</div>

</body>
</html>