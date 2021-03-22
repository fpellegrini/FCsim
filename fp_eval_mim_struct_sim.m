function fp_eval_mim_struct_sim(ip)

fp_addpath

DIRLOG ='/home/bbci/data/haufe/Franziska/log/mim_sim3/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

rng('shuffle')

%%
%prevent array jobs to start at exactly the same time
iit = str2num(getenv('SGE_TASK_ID'))

% for ip = varyParam
clear nInteractions nRegionInts SNR noise_mix nlag filtertype hemisym

[nInteractions,nRegionInts,SNR,noise_mix,nlag,filtertype] = fp_get_params(ip);

for iInt = nInteractions
    for iReg = nRegionInts
        for isnr = SNR
            for iss = noise_mix
                for ilag = nlag
                    for ifi = 1:numel(filtertype)
                        
                        ifilt = filtertype{ifi};
                        
                        iit
                        %create logfile for parallization
                        logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_iter%d'...
                            ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,iit);
                        %
                        if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
                            eval(sprintf('!touch %s%s_work',DIRLOG,logname))
                            fprintf('Working on %s. \n',logname)
                            
                            
                            params.iInt = iInt;
                            params.iReg = iReg;
                            params.isnr = isnr;
                            params.iss = iss;
                            params.ilag = ilag;
                            params.ifilt = ifilt;
                            params.iit = iit;
                            params.ip = ip;
                            params.logname = logname;
                            
                            
                            fp_data2mim_sim(params)
                            
                            eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
                            
                        end  %work_done
                        
                    end %filtertype
                end %lag
            end %noise_mix
        end %snr
    end %iReg
end %iInt

% end %varyParams


