function fp_eval_mim_struct_sim

DIRLOG = './log/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

rng('shuffle')
varyParam = 1:8;
nit = 100;

%%


for ip = varyParam
    clear nInteractions nRegionInts SNR
    
    [nInteractions,nRegionInts,SNR,noise_mix,nlag,filtertype,hemisym] = fp_get_params(ip);
    
    for iInt = nInteractions
        for iReg = nRegionInts
            for isnr = SNR
                for iss = noise_mix
                    for ilag = nlag
                        for ifilt = filtertype
                            for ihemi = hemisym
                                for iit = 1: nit
                                    
                                    
                                    %create logfile for parallization
                                    logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d_iter%d'...
                                        ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,ihemi,iit);
                                    
                                    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
                                        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
                                        fprintf('Working on %s. \n',logname)
                                        
                                        params.iInt = iInt;
                                        params.iReg = iReg;
                                        params.isnr = isnr;
                                        params.iss = iss;
                                        params.ilag = ilag;
                                        params.ifilt = ifilt;
                                        params.ihemi = ihemi;
                                        params.iit = iit; 
                                        params.ip = ip;
                                        
                                                                           
                                        fp_mim_struct_sim(params)
                                        
                                        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
                                        
                                    end  %work_done
                                end %iit
                                
                            end %hemispheres
                        end %filtertype
                    end %lag
                end %noise_mix
            end %snr
        end %iReg
    end %iInt
    
end %varyParams


