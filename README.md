# FCsim
Matlab code for the paper "Identifying best practices for detecting inter-regional functional connectivity from EEG" by Franziska Pellegrini, Arnaud Delorme, Vadim Nikulin, and Stefan Haufe. 

This is a simulation that generates artificial EEG sensor time series and compares different pipelines that aim at 
estimating source-level inter-regional functional connectivity from EEG sensor data.

The simulation is started by [fp_eval_mim_struct_sim.m](fp_eval_mim_struct_sim.m) in Matlab. The simulation on the comparison between different interaction delays starts with [fp_eval_lagsim.m](fp_eval_lagsim.m). 

The following Matlab toolboxes are required:    
- [MVGC](https://github.com/lcbarnett/ssgc) (see also [here](https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html) for an up-to-date version)       
- [Dugh-NeurIPS-2021](https://github.com/AliHashemi-ai/Dugh-NeurIPS-2021) (only for Champagne source localization)    

The recommended methods and pipelines introduced in this project have been implemented in the ROIconn plugin to the free and open source [EEGlab toolbox](https://github.com/arnodelorme/roiconnect).

The authors would be grateful if published reports of research using this code (or a modified version, maintaining a significant portion of the original code) would cite the following article: 
> Pellegrini, F., Delorme, A., Nikulin, V. & Haufe, S., 2022. Identifying best practices for detecting inter-regional functional connectivity from EEG. bioRxiv 2022.10.05.510753. https://doi.org/10.1101/2022.10.05.510753
