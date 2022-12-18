# cptp-analysis
Code for analysis of MD simulations of CPTP described in "Ceramide-1-phosphate transfer protein enhances lipid transport by disrupting hydrophobic lipid-membrane contacts" (bioRxiv DOI: [https://doi.org/10.1101/2022.09.10.507427](https://doi.org/10.1101/2022.09.10.507427)) by Jula R. Rogers and Phillip L. Geissler.

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

### Number of hydrophobic lipid–protein contacts
Property: `mem_prot_mindcc_ncc`

Number of hydrophobic contacts between protein hydrophobic carbons and lipid hydrophobic carbons, $n_{\rm CC}^{(\rm lip-prot)}$, and minimum distance between these pairs of hydrophobic carbons, $\min(d_{\rm CC}^{(\rm lip-prot)})$. Code outputs $n_{\rm CC}^{(\rm lip-prot)}$ and $\min(d_{\rm CC}^{(\rm lip-prot)})$ for each membrane lipid found in the leaflet that CPTP is bound for every frame of the trajectory and also for each protein residue of CPTP.

### Average lipid acyl chain order parameter
Property: `mem_scc`

Calculates average order parameter of each lipid's acyl chains, $S_{\rm CC}$, using [LiPyphilic](https://lipyphilic.readthedocs.io/en/latest/index.html).

## Calculation of average membrane properties in $xy$ plane
In general, average properties of lipids as a function of their displacement from CPTP in the $xy$ plane can be calculated using:

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

### Average lipid acyl chain order parameter in $xy$ plane
Property: `scc`

Average $S_{\rm CC}$ as a function of $x$ and $y$ is calculated from the output of `calc_mem_scc_allatom.py`.

### Average area per lipid in $xy$ plane
Property: `aplfatslim`

Average area per lipid, $A_{\rm lip}$, calculated with [FATSLiM](https://pythonhosted.org/fatslim/index.html), as a function of $x$ and $y$. Unlike the above properties, `<property_per_lipid_txt_file>` is not provided as input. Instead, `<fatslim_basename>` is supplied as the last command line argument. Files output from FATSLiM with $A_{\rm lip}$ for each lipid must be named `<fatslim_basename>_frame_%05d.csv`.

```
python calc_mem_av_aplfatslim_xy.py <gro> <xtc> <fatslim_basename>
```

### Average membrane thickness in $xy$ plane
Property: `thickness`

Average membrane thickness defined as the average distance in $z$ between phosphorous atoms of lipids in the top and bottom leaflets. Unlike the above properties, `<property_per_lipid_txt_file>` is not provided as input.

```
python calc_mem_av_thickness_xy.py <gro> <xtc>
```

## Calculation of membrane property distributions as a function of radial distance from CPTP

In general, distributions (average, standard deviation, and statistics for box-and-whisker plots) of each membrane property as a function of radial distance from CPTP can be calculated using:

```
python calc_mem_boxwhisker_<property>_radial.py <gro> <xtc> <property_per_lipid_txt_file>
```
where `<property_per_lipid_txt_file>` is output from analysis described above, and <gro> is a gromacs structure file and <xtc> is a gromacs trajectory file for the simulation used for the analysis.

### Lipid tail orientation
Property: `costhetaz`

Distribution of $\cos(\theta_z)$ as a function of radial distance is calculated from the output of `calc_mem_chain_orientation_allatom.py`.

### Membrane hydrophobicity
Property: `mindcc_ncc`

Distributions of $n_{\rm CC}$ and $\min(d_{\rm CC})$ as a function of radial distance are calculated from the output of `calc_mem_mindcc_ncc_allatom.py`.

### Average lipid acyl chain order parameter
Property: `scc`

Distribution of $S_{\rm CC}$ as a function of radial distance is calculated from the output of `calc_mem_scc_allatom.py`.

### Membrane thickness
Property `thickness`

Distribution of membrane thickness as a function of radial distance is calculated. Unlike the above properties, `<property_per_lipid_txt_file>` is not provided as input.

```
python calc_mem_boxwhisker_thickness_radial.py <gro> <xtc>
```

### Area per lipid
Property: `aplfatslim`

Distribution of $A_{\rm lip}$, calculated with [FATSLiM](https://pythonhosted.org/fatslim/index.html), as a function of radial distance is calculated. Instead, `<fatslim_basename>` is supplied as the last command line argument. Files output from FATSLiM with $A_{\rm lip}$ for each lipid must be named `<fatslim_basename>_frame_%05d.csv`.

```
python calc_mem_boxwhisker_aplfatslim_radial.py <gro> <xtc> <fatslim_basename>
```

## Calculation of average membrane properties within 5 Angstroms of each residue of CPTP
In general, averages (and standard deviations) of each lipid property for lipids within 5 Angstrom of a given residue of CPTP can be calculated using:

```
python calc_mem_av_<property>_protres.py  <gro> <xtc> <property_per_lipid_txt_file>
```
where `<property_per_lipid_txt_file>` is output from analysis described above, and <gro> is a gromacs structure file and <xtc> is a gromacs trajectory file for the simulation used for the analysis.

### Lipid tail orientation
Property: `costhetaz`

Average $\cos(\theta_z)$ of lipids within 5 Angstrom of each protein residue is calculated from the output of `calc_mem_chain_orientation_allatom.py`.

### Membrane hydrophobicity
Property: `mindcc_ncc`

Average $n_{\rm CC}$ and $\min(d_{\rm CC})$ of lipids within 5 Angstrom of each protein residue are calculated from the output of `calc_mem_mindcc_ncc_allatom.py`.

### Average lipid acyl chain order parameter
Property: `scc`

Average $S_{\rm CC}$ of lipids within 5 Angstrom of each protein residue is calculated from the output of `calc_mem_scc_allatom.py`.

### Area per lipid
Property: `aplfatslim`

Average $A_{\rm lip}$, calculated with [FATSLiM](https://pythonhosted.org/fatslim/index.html), of lipids within 5 Angstrom of each protein residue is calculated. Instead, `<fatslim_basename>` is supplied as the last command line argument. Files output from FATSLiM with $A_{\rm lip}$ for each lipid must be named `<fatslim_basename>_frame_%05d.csv`.

```
python calc_mem_av_aplfatslim_protres.py <gro> <xtc> <fatslim_basename>
```

## Requirements
With the exeption of the code using [LiPyphilic](https://lipyphilic.readthedocs.io/en/latest/index.html), all code requires:

* Python 2.7+
* Numpy 1.16+
* MDAnalysis 0.20

Code using [LiPyphilic](https://lipyphilic.readthedocs.io/en/latest/index.html), specifically `calc_mem_scc_allatom.py` requires:

* Python 3.7+
* Numpy 1.21+
* MDAnalysis 2.1
