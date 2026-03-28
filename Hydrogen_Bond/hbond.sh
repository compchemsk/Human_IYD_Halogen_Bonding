for np in 3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 9.8 10.8 11.8 12.8 13.8 14.8 15.8 16.8 17.8 18.8
do
cd US-${np}
rm -rf hbond_percentages-${np}.txt
echo "
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from collections import Counter

u = mda.Universe('4ttc-subs-EP.parm7', '${np}-us-10ns.nc')

hbond_analysis = HydrogenBondAnalysis(
    universe=u,
    donors_sel='resid 1-442 and (name N* or name O*)',
    hydrogens_sel='resid 1-442 and name H*',
    acceptors_sel='resid 1-442 and (name O* or name N*)',
    d_h_cutoff=1.2,
    d_a_cutoff=3.0,
    d_h_a_angle_cutoff=150,
    update_selections=True
)

hbond_analysis.run()

hbonds_data = hbond_analysis.results.hbonds
hbond_counts = Counter()
seen_pairs_per_frame = dict()

for hbond in hbonds_data:
    frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle = hbond
    donor_idx = int(donor_idx)
    acceptor_idx = int(acceptor_idx)

    donor_atom = u.atoms[donor_idx]
    acceptor_atom = u.atoms[acceptor_idx]

    donor_resid = donor_atom.resid
    acceptor_resid = acceptor_atom.resid

    if (donor_resid == 442 and 1 <= acceptor_resid <= 441) or (acceptor_resid == 442 and 1 <= donor_resid <= 441):
        partner_resid = donor_resid if acceptor_resid == 442 else acceptor_resid
        key = (int(frame), partner_resid)

        if key not in seen_pairs_per_frame:
            seen_pairs_per_frame[key] = True
            hbond_counts[partner_resid] += 1

total_frames = hbond_analysis.n_frames
percentage_hbonds = {resid: (count / total_frames) * 100 for resid, count in hbond_counts.items()}

with open(f'hbond_percentages-${np}.txt', 'w') as out:
    out.write('Residue_ID\\tHbond_Percentage\\n')
    for resid in sorted(percentage_hbonds):
        out.write(f'{resid}\\t{percentage_hbonds[resid]:.2f}\\n')

print(f'Hydrogen bond analysis complete for distance ${np}. Results saved to hbond_percentages-${np}.txt')
" > ${np}.py

python ${np}.py
rm ${np}.py
sed -i '1d' hbond_percentages-${np}.txt
cat hbond_percentages-${np}.txt | sort -nrk2 > ../hbond-${np}.dat
cd ../
sleep 1s
done
paste $(for np in 3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 9.8 10.8 11.8 12.8 13.8 14.8 15.8 16.8 17.8 18.8; do echo hbond-${np}.dat; done) > hbond_all.dat
rm hbond-*dat

