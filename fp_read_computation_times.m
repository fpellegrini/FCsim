function fp_read_computation_times

fp_addpath

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim5/';
DIROUT_ctimes = '/home/bbci/data/haufe/Franziska/data/mim_sim5/computation_times/';
if ~exist(DIROUT_ctimes); mkdir(DIROUT_ctimes); end

%%
for ip = 1
    
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
                
                                try
                                    ta(1,iit) = t.signal;
                                catch
                                    ta(1,iit) = nan;
                                end
                                try
                                    ta(2,iit) = t.snr;
                                catch
                                    ta(2,iit) = nan;
                                end
                                try
                                    ta(3,iit) = t.filter;
                                catch
                                    ta(3,iit) = nan;
                                end
                                try
                                    ta(4:numel(t.pips)+3,iit) = t.pips;
                                catch
                                    ta(4,iit) = nan;
                                end
                                
                                try
                                    corr_voxmim_all(iit,:) = corr_voxmim;
                                    corr_voxmic_all(iit,:) = corr_voxmic;
                                end
                                
                                for oo = params.pips
                                    try
                                        varex_all(iit,oo) = mean(to_save{oo}.varex);
                                    end
                                end
                                
                            end  %iit
                            
                            info = {'signal','snr','filter','pips'};
                            pips = params.pips; 
                            
                            outname = [DIROUT_ctimes logname '.mat'];
                            save(outname,'ta','info','pips','corr_voxmim_all','corr_voxmic_all', 'varex_all','-v7.3')
                            
                        end %filtertype
                    end %lag
                end %noise_mix
            end %snr
        end %iReg
    end %iInt
end%ip


%%

% figure; 
% figone(10,12)
% b=bar(varex(1:6));
% grid on 
% b.FaceColor = [0.8 0.7 0.6];
% xlabel('Number of PCs')
% ylabel('Variance explained')