function fp_read_computation_times

fp_addpath

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim5/';
DIROUT_ctimes = '/home/bbci/data/haufe/Franziska/data/mim_sim5/computation_times/';
if ~exist(DIROUT_ctimes); mkdir(DIROUT_ctimes); end

%%
for ip = [2:5 7]
    
    clear nInteractions nRegionInts SNR noise_mix nlag filtertype hemisym    
    [nInteractions,nRegionInts,SNR,noise_mix,nlag,filtertype, dimred] = fp_get_params(ip);
    
    for iInt = nInteractions
        for iReg = nRegionInts
            for isnr = SNR
                for iss = noise_mix
                    for ilag = nlag
                        for ifi = 1:numel(filtertype)
                            
                            clear ifilt
                            ifilt = filtertype{ifi};
                            
                            for iit = 1:100
                                iit
                                %create logfile for parallelization
                                logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_iter%d'...
                                    ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,dimred,iit);
                                
                                load([DIRIN 'mim_' logname '.mat'])
        
                                ta(1,iit) = t.signal;
                                ta(2,iit) = t.snr;
                                ta(3,iit) = t.filter;
                                ta(4:numel(t.pips)+3,iit) = t.pips;      
                                
                            end  %iit
                            
                            info = {'signal','snr','filter','pips'};
                            pips = params.pips; 
                            
                            outname = [DIROUT_ctimes logname '.mat'];
                            save(outname,'ta','info','pips', '-v7.3')
                            
                        end %filtertype
                    end %lag
                end %noise_mix
            end %snr
        end %iReg
    end %iInt
end%ip