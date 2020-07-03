function [nInteractions, nRegionInts,SNR,noise_mix,nlag,filtertype,hemisym] = fp_get_params(ip)

if ip == 1
    %defaults
    nInteractions = 1;
    nRegionInts = 2;
    SNR = 0.5;
    noise_mix = 0.5;
    nlag = 1;
    filtertype= 'd'; %dics
    hemisym = 0;
    
elseif ip == 2
    %vary nInteractions
    nInteractions = 2:5;
    nRegionInts = 1;
    SNR = 0.5;
    noise_mix = 0.5;
    nlag = 1;
    filtertype= 'd'; %dics
    hemisym = 0;
    
elseif ip == 3
    %vary nRegionInts
    nInteractions = 1;
    nRegionInts = 2:3;
    SNR = 0.5;
    noise_mix = 0.5;
    nlag = 1;
    filtertype= 'd'; %dics
    hemisym = 0;
    
elseif ip == 4
    %vary SNR
    nInteractions = 1;
    nRegionInts = 1;
    SNR = [0.1:0.1:0.4, 0.6:0.1:0.9];
    noise_mix = 0.5;
    nlag = 1;
    filtertype= 'd'; %dics
    hemisym = 0;
    
elseif ip == 5
    %vary signal strength
    nInteractions = 1;
    nRegionInts = 1;
    SNR = 0.5;
    noise_mix = [0.1:0.1:0.4 0.6:0.1:0.9];
    nlag = 1;
    filtertype= 'd'; %dics
    hemisym = 0;
    
elseif ip == 6
    %vary lag size
    nInteractions = 1;
    nRegionInts = 1;
    SNR = 0.5;
    noise_mix = 0.5;
    nlag = 2; %small (0 to 5 samples (=1)) or large (5 to 20 samples (=2))
    filtertype= 'd'; %dics
    hemisym = 0;
    
elseif ip == 7
    %vary filter
    nInteractions = 1;
    nRegionInts = 1;
    SNR = 0.5;
    noise_mix = 0.5;
    nlag = 1;
    filtertype= 'e'; %eloreta
    hemisym = 0;
    
elseif ip == 8
    %vary hemisphere symmetry
    nInteractions = 1;
    nRegionInts = 1;
    SNR = 0.5;
    noise_mix = 0.5;
    nlag = 1;
    filtertype= 'd'; %dics
    hemisym = 1; %symmetrize hemispheres
end