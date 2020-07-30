function mim_struct_sim

fp_addpath

DIROUT = '/home/bbci/data/haufe/Franziska/data/';
DIRLOG = './log/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

rng('shuffle')
varyParam = 1:8;
nit = 1; %500;
fres = 40;
n_trials = 200;

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
                                
                                %in each iteration, a new signal with new
                                %interacting sources and voxels is generated
                                for iit = 1: nit
                                    
                                    %create logfile for parallization
                                    logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d_iter%d'...
                                        ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,ihemi,iit);
                                    
                                    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
                                        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
                                        fprintf('Working on %s. \n',logname)
                                        


                                        
                                        %% signal generation
                                        
                                        clearvars -except DIROUT DIRLOG ...
                                            iit nit nInteractions nRegionInts SNR noise_mix ...
                                            nlag filtertype hemisym ip varyParam...
                                            fres n_trials iss iInt iReg isnr ilag ifilt ihemi
                                        
                                        params.iInt = iInt;
                                        params.iReg = iReg;
                                        params.isnr = isnr;
                                        params.iss = iss;
                                        params.ilag = ilag;
                                        params.ifilt = ifilt;
                                        params.ihemi = ihemi;
                                        
                                        
                                        % ROI labels
                                        % In D.sub_ind_roi, there are the randomly
                                        % selected voxels of each region
                                        fprintf('Getting atlas positions... \n')
                                        tic
                                        D = fp_get_Desikan(params.iReg);
                                        toc
                                        
                                        %signal generation
                                        fprintf('Signal generation... \n')
                                        tic
                                        [signal_sensor,gt,L,iroi_seed, iroi_tar] = fp_generate_mim_signal(params, ...
                                            fres,n_trials, D);
                                        toc
                                        
                                        
                                        %% get CS and filter A
                                        %parameters
                                        id_trials_1 = 1:n_trials;
                                        id_trials_2 = 1:n_trials;
                                        id_meg_chan = 1:size(signal_sensor,1);
                                        nmeg = numel(id_meg_chan);
                                        regu=.000001;
                                        
                                        %cross spectrum
                                        fprintf('Calculating cross spectrum... \n')
                                        tic
                                        CS = fp_tsdata_to_cpsd(signal_sensor,fres,'WELCH',...
                                            [id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
                                        CS(:,:,1)=[];
                                        nfreq = size(CS,3);
                                        toc
                                        
                                        %leadfield
                                        L3 = L(:, D.ind_cortex, :);
                                        for is=1:D.nvox
                                            clear L2
                                            L2 = L3(:,is,:);
                                            
                                            %remove radial orientation
                                            clear u s
                                            [u, s, v] = svd(squeeze(L2));
                                            L_backward(:,is,:) = u(:,:)*s(:,1:2);
                                        end
                                        ni = size(L_backward,3);
                                        
                                        %construct source filter
                                        if strcmp(ifilt,'e')
                                            A = squeeze(mkfilt_eloreta_v2(L_backward));
                                            A = permute(A,[1, 3, 2]);
                                            fqA = ones(1,nfreq);%only one filter for all freqs.
                                            nfqA = 1;
                                            
                                        elseif strcmp(ifilt,'d')
                                            
                                            A=zeros(nmeg,ni,D.nvox,nfreq);
                                            
                                            for ifrq = 1:nfreq
                                                cCS = CS(:,:,ifrq);
                                                lambda = mean(diag(real(cCS)))/100;
                                                
                                                CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
                                                
                                                for is=1:D.nvox %iterate across nodes
                                                    Lloc=squeeze(L_backward(:,is,:));
                                                    A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
                                                end
                                            end
                                            fqA = 1:nfreq; %This filter is frequency specific.
                                            nfqA = nfreq;
                                            
                                            
                                        elseif strcmp(filtertype,'l')                 
                                        end
                                        
                                        %% calculate MIM
                                        %pca pipeline ('all' 8 pipelines + baseline)
                                        [mic, mim] = fp_get_mim(A,CS,fqA,D,params.ihemi,'all');
                                        
                                        %% performance measures
                                        fprintf('Performance measures... \n') 
                                        tic
                                        [PERFORMANCE, BASELINE] = fp_get_performance(gt, mic, mim, params); 
                                        toc
                                        
                                        %% save performance and baseline 
                                        fprintf('Saving... \n')
                                        tic
                                        outname = sprintf('%smim_%s.mat',DIROUT,logname);
                                        save(outname,'PERFORMANCE','BASELINE','-v7.3')
                                        toc
                                        
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
