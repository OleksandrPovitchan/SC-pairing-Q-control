# SC-pairing-Q-control
***(numerical analysis for arXiv:2503.20426)***

One can find the list of packages in the conda environment used when executing Python scripts in the file *package_list.txt*

**Brief description of the provided Python files**

All the .py files starting with 'grid' produce data for plotting *FIG. 2.* from the article. The process for the uncontrolled evolution (UE) is described below. Similar actions are required to simulate data in the case of Lyapunov and asymptotic controls (LQC & AQC). Suppose one wants to calculate the final and maximal values of the order parameter over the region from 0 to 60 in the angular frequency and from 0 to 3 in the amplitude, and wants to split the process into 12 parts. In that case, one can create one of 12 smaller grids by setting the grid parameters in *parameters.py* and running the following sequence of scripts: *grid-UE.py* [or, equivalently, *grid-UE-to-state.py* -> *grid-UE-from-state.py* (this approach requires a lot of storage, since files containing the whole evolution of the wave function might weight several GB already for 8 sites)] -> *grid-UE-processing.py*. When the process is completed for all the smaller grids, one should run *grid-UE-global-dict.py* and proceed with files in the directory *analysis*.

*UE-eta_sq.py* produces data for *UE-fig-1-b.py* (the latter creates *FIG. 1. b.*). *UE-eta_sq.py* + *UE-spectral-weights.py* together produce data for *UE-fig-1-a.py* (the latter creates *FIG. 1. a.*)

*FIG_3-UE-LC.py* + *FIG_3-UE-AC.py* together produce data for *FIG_3-a.py* (the latter creates *FIG. 3. a.*). *FIG_3-UE-LC.py* + *FIG_3-spectral-weights.py* together produce data for *FIG_3-b.py* (the latter creates *FIG. 3. b.*)

*FIG_3-UE-LC.py* + *FIG_4-LC-saturation.py* together produce data for *FIG_4-a.py* (the latter creates *FIG. 4. a.*). By adding to the previous two files *FIG_4-19.1-0.2-asympt-sat.py* and *FIG_4-17.0-0.4-asympt-sat.py*, we prepare data for *FIG_4-b.py* (the latter creates *FIG. 4. b.*).
