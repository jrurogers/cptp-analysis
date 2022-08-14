# cptp-analysis
Code for analysis of MD simulations of CPTP described in "Ceramide-1-phosphate transfer protein enhances lipid transport by disrupting hydrophobic lipid-membrane contacts"

## Analysis of simulation trajectories
General usage of analysis code:
```
python calc_<property>_<simulation_type>.py <gro> <xtc>
```
where `<property>` is the property to calculate (see below descriptions), `<simulation_type>` indicates the resolution of the simulation to analyze (`allatom` for all-atom simulations using the CHARMM36 force field or `martini` for coarse-grained simulations using the MARTINI2.3P force field), `<gro>` is a gromacs structure file for the simulation, and `<xtc>` is a gromacs trajectory file for the simulation. Analysis that requires additional inputs is described below.

### Insertion depth
Property: `insertion_depth`

The insertion depth of each residue is the average signed distance in $z$, which is the axis normal to the membrane, between the center of mass (COM) of each residue and the average position of phosphate atoms (beads) in all-atom (coarse-grained) simulations.

### Orientation of helix $\alpha2$
Property: `angle_phi_alpha6plane`

The orientation of helix $\alpha2$ is quantified by its polar angle $\theta_{\alpha 2}$ and azimuthal angle $\phi_{\alpha 2}$ in a coordinate system in which the $xy$ plane is approximately the membrane surface. Code outputs $\theta_{\alpha 2}$ and $\phi_{\alpha 2}$ for every frame in the trajectory in addition to average values and 1D & 2D probability distributions.

### Lipid tail orientation
Property: `mem_chain_orientation`

Orientations of the tails of each lipid in the membrane are quantified by their angle $\theta_z$ with the $z$-axis. Code outputs $\theta_z$ of both tails of each membrane lipid found in the leaflet that CPTP is bound for each frame in the trajectory.

### Membrane hydrophobicity
Property: `mem_mindcc_ncc`

Membrane hydrophobicity quantified by the number of hydrophobic contacts a lipid makes with surrounding membrane lipids, $n_{\rm CC}$, and minimum distance between a pair of hydrophobic carbons on two different lipids, $\min(d_{\rm CC})$. Code outputs $n_{\rm CC}$ and $\min(d_{\rm CC})$ for each membrane lipid found in the leaflet that CPTP is bound for every frame of the trajectory.

### RMSD to reference structure
Property: `rmsd2pdb`

Root-mean-square deviation of ${\rm C}\alpha$ calculated between simulated configurations of CPTP and a reference PDB structure provided as last command line argument:
```
python calc_rmsd2pdb_<simulation_type>.py <gro> <xtc> <pdb>
```
PDBs used for analysis provided in `reference_structures/`: For all-atom simulations, `4k85_ca.pdb`, which includes only ${\rm C}\alpha$ atoms from the crystal structure provided in PDB 4K85, was used. For coarse-grained simulations, `4k85_bb_martini.pdb`, which includes only ${\rm C}\alpha$ atoms from the crystal structure provided in PDB 4K85 renamed to BB, was used. Code outputs RMSD for every frame of the trajectory.

### RMSD between C1P bound and apo forms
Property: `rmsd_apo2c1p`

Root-mean-square deviation of ${\rm C}\alpha$ atoms between simulated C1P bound and apo forms of CPTP:

$$ \text{RMSD}\_{{\rm C}\alpha}({\rm C1P} - {\rm apo}) = \sqrt{t^{-2} \sum_{ij} (\vec{r}\_{\rm C1P}(t_i) - \vec{r}\_{\rm apo}(t_j))^2} $$

where $i$ and $j$ index all $t$ frames in each trajectory and $\vec{r}\_{\rm C1P}(t)$ and $\vec{r}\_{\rm apo}(t)$ are the ${\rm C}\alpha$'s positions in the C1P-bound and apo forms, respectively. Prior to calculation of $\text{RMSD}\_{{\rm C}\alpha}({\rm C1P} - {\rm apo})$, the simulation trajectories are both aligned to a common reference frame. The reference structure used for analysis is provided in `reference_structures/`: `rmsd_apo2c1p_ref.gro`, which is the first frame from the all-atom solution simulations of the C1P bound form of CPTP. This reference structure, gromacs structure and trajectories files for simulations of the apo and C1P bound forms are provided as command line arguments:
```
python calc_rmsd_apo2c1p_allatom.py <ref_gro> <apo_gro> <apo_xtc> <c1p_gro> <c1p_xtc>
```

### RMSF
Property: `rmsf`

Root-mean-square fluctuation calculated for ${\rm C}\alpha$ atoms.


## Calculation of average membrane properties in $xy$ plane
Average properties of lipids as a function of their displacement from CPTP in the $xy$ plane can be calculated using:
```
python calc_mem_av_<property>_xy.py <gro> <xtc> <property_per_lipid_txt_file>
```
where `<property_per_lipid_txt_file>` is output from analysis described above, and `<gro>` is a gromacs structure file and `<xtc>` is a gromacs trajectory file for the simulation used for the analysis. A reference structure (`reference_structures/av_xy_plane_ref.gro`) is used for alignment and centering of CPTP's COM at the box center prior to calculating each lipid's $xy$ position relative to CPTP.

### Average lipid orientation in $xy$ plane
Property: `costhetaz`

Average $\cos(\theta_z)$ as a function of $x$ and $y$ is calculated from the output of `calc_mem_chain_orientation_allatom.py`.

### Average membrane hydrophobicity in $xy$ plane
Property: `mindcc_ncc`

Average $n_{\rm CC}$, $\min(d_{\rm CC})$, and $r_{\rm LxS}$ as a function of $x$ and $y$ are calculated from the output of `calc_mem_mindcc_ncc_allatom.py`. The reaction coordinate for passive lipid transport via solvent, $r_{\rm LxS}$, is a function of $n_{\rm CC}$ and $\min(d_{\rm CC})$ as described in ["Membrane hydrophobicity determines the activation free energy of passive lipid transport"](https://www.sciencedirect.com/science/article/abs/pii/S0006349521006020) and ["Breakage of hydrophobic contacts limits the rate of passive lipid exchange between membranes"](https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.0c04139). The coefficients used to calculate $r_{\rm LxS}$ are set in the code for all-atom simulations and would need to be changed for analysis of coarse-grained simulations.

## Requirements

* Python 2.7+
* Numpy 1.16+
* MDAnalysis 0.20
