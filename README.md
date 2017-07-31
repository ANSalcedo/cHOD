## cHOD
cHOD takes an input hdf5 halo catalog, as well as environmental densities and outputs a hdf5 galaxy catalog (positions, velocities, sat/cen flag, and host ID) using a traditional HOD framework (with velocity bias) extended to include environment-dependent galaxy assembly bias in addition to allowing the satellite profile to deviate from the host's NFW profile.

## Usage
./compute_mocks Omega_m0 siglogM logMmin logM0 logM1 alpha q_env del_gamma alpha_cen alpha_sat boxsize seed halo_catalog_file galaxy_mock_file halo_environment_file

## HOD Prescription
We use a traditional HOD prescription:

![Alt Text](https://github.com/ANSalcedo/Emulator-Pipeline/raw/master/cHOD/readmegifs/standard_HOD.gif)

while also allowing for environment-dependent galaxy assembly bias. This is implemented by allowing the expectation of hosting a central to depend, at fixed mass, on overdensity percentile:

![Alt Text](https://github.com/ANSalcedo/Emulator-Pipeline/raw/master/cHOD/readmegifs/Q_param.gif)

We also allow the profile of the galaxies to deviate from the NFW profile of the host halo:

![Alt Text](https://github.com/ANSalcedo/Emulator-Pipeline/raw/master/cHOD/readmegifs/delgam.gif)

We implement the velocity bias model of Guo et al. 2015 (http://adsabs.harvard.edu/abs/2015MNRAS.453.4368G) where the central velocity is given by:

![Alt Text](https://github.com/ANSalcedo/Emulator-Pipeline/raw/master/cHOD/readmegifs/alpha_cen.gif)

and the satellite velocities are given by:

![Alt Text](https://github.com/ANSalcedo/Emulator-Pipeline/raw/master/cHOD/readmegifs/alpha_sat.gif)
