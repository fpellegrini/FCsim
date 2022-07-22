function fp_eval_mim_struct_sim(ip)

fp_addpath

%logfile to faciliatate parallel computing on the cluster 
DIRLOG ='/home/bbci/data/haufe/Franziska/log/mim_sim5/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

rng('shuffle')

%%
%iteration number equals cluster job number 
iit = str2num(getenv('SGE_TASK_ID'));

clear nInteractions nRegionInts SNR noise_mix nlag filtertype hemisym

%get parameters for experimental setup 
[nInteractions,nRegionInts,SNR,noise_mix,nlag,filtertype, dimred] = fp_get_params(ip); 

for iInt = nInteractions
    for iReg = nRegionInts
        for isnr = SNR
            for iss = noise_mix
                for ilag = nlag
                    for ifi = 1:numel(filtertype)
                        
                        clear ifilt
                        ifilt = filtertype{ifi};
                        
                        iit
                        %create logname
                        if ip==9
                            logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_corr_iter%d'...
                                ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,dimred,iit);
                        else
                            logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_%s_iter%d'...
                                ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,dimred,iit);
                        end

                        if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
                            eval(sprintf('!touch %s%s_work',DIRLOG,logname)) %create logfile for parallelization
                            fprintf('Working on %s. \n',logname)                           
                            
                            %setup params struct with specificaltions for
                            %experimental setup 
                            params.iInt = iInt; %number of interactions 
                            params.iReg = iReg; %number of active voxels per region 
                            params.isnr = isnr; %signal-to-noise ratio
                            params.iss = iss; %noise mix
                            params.ilag = ilag; %magitude of time delay (ilag=2 -> time delay between 50 and 200 Hz)
                            params.ifilt = ifilt; %source projection filter type 
                            params.dimred = dimred; %type of dimensionality reduction (pca vs ssd, but ssd not properly investigated here) 
                            params.iit = iit;
                            params.ip = ip;
                            params.logname = logname;
                            
                            if strcmp(params.ifilt,'d') 
                                %dics pipeline (in frequency space from the beginning)
                                fp_data2mim_sim_dics(params) 
                            else
                                fp_data2mim_sim(params)
                            end
                            
                            eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
                            
                        end  %work_done
                        
                    end %filtertype
                end %lag
            end %noise_mix
        end %snr
    end %iReg
end %iInt



