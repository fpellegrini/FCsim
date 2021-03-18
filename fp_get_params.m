function [nInteractions, nRegionInts,SNR,noise_mix,nlag,filtertype] = fp_get_params(ip)

if ip == 1
    %defaults
    nInteractions = 2;
    nRegionInts = 1;
    SNR = 0.7;
    noise_mix = 0.5;
    nlag = 2;
    filtertype= {'l'}; %lcmv
    
elseif ip == 2
    %vary nInteractions
    nInteractions = [1 3:5];
    nRegionInts = 1;
    SNR = 0.7;
    noise_mix = 0.5;
    nlag = 2;
    filtertype= {'l'}; %lcmv
    
elseif ip == 3
    %vary nRegionInts
    nInteractions = 2;
    nRegionInts = 2;
    SNR = 0.7;
    noise_mix = 0.5;
    nlag = 2;
    filtertype= {'l'}; %lcmv
    
elseif ip == 4
    %vary SNR
    nInteractions = 2;
    nRegionInts = 1;
    SNR = [0.1 0.3 0.5 0.9];
    noise_mix = 0.5;
    nlag = 2;
    filtertype= {'l'}; %lcmv
    
elseif ip == 5
    %vary noise_mix
    nInteractions = 2;
    nRegionInts = 1;
    SNR = 0.7;
    noise_mix = [0 0.25 0.75 1];
    nlag = 2;
    filtertype= {'l'}; %lcmv
    
elseif ip == 6
    %vary lag size
    nInteractions = 2;
    nRegionInts = 1;
    SNR = 0.7;
    noise_mix = 0.5;
    nlag = 1; %small (0 to 5 samples (=1)) or large (5 to 20 samples (=2))
    filtertype= {'l'}; %lcmv
    
elseif ip == 7
    %vary filter
    nInteractions = 2;
    nRegionInts = 1;
    SNR = 0.7;
    noise_mix = 0.5;
    nlag = 2;
    filtertype= {'e','c'};
    
end