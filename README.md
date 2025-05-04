# SC-pairing-Q-control
***(numerical analysis for arXiv:2503.20426)***

**Brief description of the provided Python files**

All the .py files starting with 'grid' produce data for plotting *FIG. 2.* from the article. The process for the uncontrolled evolution (UE) is described below. Similar actions are required to simulate data in the case of Lyapunov and asymptotic control (LQC & AQC). Suppose one wants to calculate the final and maximal values of the order parameter over the region from 0 to 60 in the angular frequency and from 0 to 3 in the amplitude, and wants to split the process into 12 parts. In that case, one can create one of 12 smaller grids by setting the grid parameters in *supporting/tools.py* and running the following sequence of scripts: *grid-UE.py* [or, equivalently, *grid-UE-to-state.py* -> *grid-UE-from-state.py* (this approach requires a lot of storage space, since files containing the whole evolution of the wave function might weight several GB already for 8 sites)] -> *grid-UE-processing.py*. When the process is completed for all the smaller grids, one should run *grid-UE-global-dict.py* and proceed with files in the directory *analysis*
