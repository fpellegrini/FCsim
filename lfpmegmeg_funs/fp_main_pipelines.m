function fp_main_pipelines 

%% pipeline lfp--meg coherence 
DIROUT = '/home/bbci/data/haufe/Franziska/data/';
% DIRLOG = '/home/bbci/data/haufe/Franziska/log/filters_dics/';

fp_savefilter([],DIROUT, DIRLOG)

fp_surrogate_coh([],DIROUT,DIRLOG)
fp_savefilter_eloreta([],DIROUT,DIRLOG)

fp_surrogate_coh_e([],DIROUT,DIRLOG)
fp_get_thresholds_ss_freq([],'imag',DIROUT)

fp_cluster_ss_c_freq([],'imag',DIROUT)
fp_cluster_g_c_freq('imag',DIROUT)

view_clusters


%% pipeline megmeg coherence 

fp_megmeg_pipeline

%% pipeline lfp meg gc 

fp_gc_pipeline1
