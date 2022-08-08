# FC_sim
Matlab code for the paper "Identifying best practices for detecting inter-regional functional connectivity from EEG" by Franziska Pellegrini, Arnaud Delorme, Vadim Nikulin, and Stefan Haufe. 

Simulation that generates artificial EEG sensor time series and compares different pipelines that aim at 
estimating source-level inter-regional functional connectivity from EEG sensor data.

The simulation is started by fp_eval_mim_struct_sim.m in Matlab. The simulation on the comparison between different interaction delays starts with fp_eval_lagsim. 

The following Matlab toolboxes are required:    
MVGC: https://github.com/lcbarnett/ssgc (see also here for an uptodate version: https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html)       
Dugh-NeurIPS-2021: https://github.com/AliHashemi-ai/Dugh-NeurIPS-2021 (only for Champagne source localization)    

The recommended methods and pipelines introduced in this have been implemented in the ROIconn plugin to the free and open source EEGlab toolbox: https://github.com/arnodelorme/roiconnect .

The authors would be grateful if published reports of research using this code (or a modified version, maintaining a significant portion of the original code) would cite the following article: tbn
