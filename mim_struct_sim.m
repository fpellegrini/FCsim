%todo: 
%1) replace tsdata with Guido's code 
%2) Performance measure 2) with whole matrix,
%2b)measure 1) auch noch über freqs aufsummiert 
%3) vary lag as a separate parameter:  random zwischen 0-5 oder 5-20
%4) Pipeline 6): rank of sum(real(CS),3) = npcs
%5) Change symmetry of hemispheres 
%6) calculate pipeline 8) region-wise 
%7) Noise defaults ändern
%8) vary eloreta 
%9) vary hemispheres symmetric or not 

%%

DIROUT = [];
DIRLOG = [];

rng('shuffle')
varyParam = 1:5;
nit = 500;
fres = 40;
n_trials = 200;

%%

for ip = varyParam
    clear nInteractions nRegionInts SNR
    
    if ip == 1
        %defaults
        nInteractions = 1;
        nRegionInts = 1;
        SNR = 0.5;
        signal_strength = 0.5;
        
    elseif ip == 2
        %vary nInteractions
        nInteractions = 2:5;
        nRegionInts = 1;
        SNR = 0.5;
        signal_strength = 0.5;
        
    elseif ip == 3
        %vary nRegionInts
        nInteractions = 1;
        nRegionInts = 2:3;
        SNR = 0.5;
        signal_strength = 0.5;
        
    elseif ip == 4
        %vary SNR
        nInteractions = 1;
        nRegionInts = 1;
        SNR = 0.1:0.1:0.9;
        signal_strength = 0.5;
        
    elseif ip == 5
        %vary SNR
        nInteractions = 1;
        nRegionInts = 1;
        SNR = 0.5;
        signal_strength = 0.1:0.1:0.9;
    end
    
    
    for iInt = nInteractions
        for iReg = nRegionInts
            for isnr = SNR
                for iss = signal_strength
                    
                    %create logfile for parallization
                    logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d',iInt,iReg,isnr*10,iss*10);
                    
                    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
                        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
                        fprintf('Working on %s. \n',logname)                       
                        
                        clear params PERFORMANCE
                        params.iInt = iInt;
                        params.iReg = iReg;
                        params.isnr = isnr;
                        params.iss = iss;
                        
                        PERFORMANCE = zeros(2,8,4,nit);
                        
                        %in each iteration, a new signal with new
                        %interacting sources and voxels is generated 
                        for iit = 1: nit
                            tic
                            %% signal generation
                            
                            clearvars -except mm_gt mc_gt bmm_gt bmc_gt GT MIM MIC ...
                                params iit nit nInteractions nRegionInts SNR ip varyParam...
                                fres n_trials PERFORMANCE iss iInt iReg isnr 
                            
                            
                            % ROI labels
                            % In D.sub_ind_roi, there are the randomly
                            % selected voxels of each region
                            D = fp_get_Desikan(params.iReg);
                            
                            %signal generation
                            [signal_sensor,gt,L,iroi_seed, iroi_tar] = fp_generate_mim_signal(params, ...
                                fres,n_trials, D);
                            
                            
                            %% get CS and filter A
                            %parameters
                            id_trials_1 = 1:n_trials;
                            id_trials_2 = 1:n_trials;
                            id_meg_chan = 1:size(signal_sensor,1);
                            nmeg = numel(id_meg_chan);
                            filtertype= 'd'; %dics
                            regu=.000001;
                            
                            %cross spectrum 
                            CS = fp_tsdata_to_cpsd(signal_sensor,fres,'WELCH',...
                                [id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
                            CS(:,:,1)=[];
                            nfreq = size(CS,3);
                            
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
                            if strcmp(filtertype,'e')
                                A = squeeze(mkfilt_eloreta_v2(L_backward));
                                A = permute(A,[1, 3, 2]);
                                fqA = ones(1,nfreq);%only one filter for all freqs.
                                nfqA = 1;
                                
                            elseif strcmp(filtertype,'d')
                                
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
                            
                            %pca pipeline ('all' 7 pipelines)
                            [mic_pca, mim_pca] = fp_get_mim_pca(A,CS,fqA,D,'all');
                            
                            %case two pipeline: mim only to reduce dimensions to 1,
                            %sum up mim and mic within regions
                            %TAKES A LOT OF MEMORY AND TIME!
                            [mic_2, mim_2, benchmark] = fp_get_mim_case2(A,CS,fqA,D);
                            
                            
                            
                            %% performance measures
                            
                            % (1)correlation mim/mic and ground truth
                            for ii = 1:5
                                cc.fixed(ii) = corr(mim_pca.fixed{ii}(:),gt(:));
                                cm.fixed(ii) = corr(mim_pca.fixed{ii}(:),gt(:));
                            end
                            cc.max = corr(mic_pca.max(:),gt(:));
                            cc.percent = corr(mic_pca.percent(:),gt(:));
                            cc.c2 = corr(mic_2(:),gt(:));
                            cm.max = corr(mim_pca.max(:),gt(:));
                            cm.percent = corr(mim_pca.percent(:),gt(:));
                            cm.c2 = corr(mim_2(:),gt(:));
                            
                            % dimensions: (mim/mic, pipeline,perfomance measure,iit)
                            PERFORMANCE(1,1:5,1,iit) = cc.fixed;
                            PERFORMANCE(2,1:5,1,iit) = cm.fixed;
                            PERFORMANCE(1,6,1,iit) = cc.max;
                            PERFORMANCE(2,6,1,iit) = cm.max;
                            PERFORMANCE(1,7,1,iit) = cc.percent;
                            PERFORMANCE(2,7,1,iit) = cm.percent;
                            PERFORMANCE(1,8,1,iit) = cc.c2;
                            PERFORMANCE(2,8,1,iit) = cm.c2;
                            
                            BENCHMARK(1,iit) = corr(benchmark(:),gt(:));
                            
                            
                            
                            % (2) correlation maxima of mim/mic and ground truth
                            gt_flat = sum(sum(gt,2),3);
                            for ii = 1:5
                                %mic
                                clear m m_max
                                m = sum(sum(mic_pca.fixed{ii},3),2);
                                m_max = fp_get_nmaxima(m,params.iInt*2);
                                ccm.fixed(ii) = corr(m_max,gt_flat);
                                
                                %mim
                                clear m m_max
                                m = sum(sum(mim_pca.fixed{ii},3),2);
                                m_max = fp_get_nmaxima(m,params.iInt*2);
                                cmm.fixed(ii) = corr(m_max,gt_flat);
                            end
                            
                            %max pipeline
                            %mic
                            clear m m_max
                            m = sum(sum(mic_pca.max,3),2);
                            m_max = fp_get_nmaxima(m,params.iInt*2);
                            ccm.max = corr(m_max,gt_flat);
                            %mim
                            clear m m_max
                            m = sum(sum(mim_pca.max,3),2);
                            m_max = fp_get_nmaxima(m,params.iInt*2);
                            cmm.max = corr(m_max,gt_flat);
                            
                            %percent pipeline
                            %mic
                            clear m m_max
                            m = sum(sum(mic_pca.percent,3),2);
                            m_max = fp_get_nmaxima(m,params.iInt*2);
                            ccm.percent = corr(m_max,gt_flat);
                            %mim
                            clear m m_max
                            m = sum(sum(mim_pca.percent,3),2);
                            m_max = fp_get_nmaxima(m,params.iInt*2);
                            cmm.percent = corr(m_max,gt_flat);
                            
                            %case2 pipeline
                            %mic
                            clear m m_max
                            m = sum(sum(mic_2,3),2);
                            m_max = fp_get_nmaxima(m,params.iInt*2);
                            ccm.c2 = corr(m_max,gt_flat);
                            %mim
                            clear m m_max
                            m = sum(sum(mim_2,3),2);
                            m_max = fp_get_nmaxima(m,params.iInt*2);
                            cmm.c2 = corr(m_max,gt_flat);
                            
                            %benchmark
                            clear m m_max
                            m = sum(sum(benchmark,3),2);
                            m_max = fp_get_nmaxima(m,params.iInt*2);
                            BENCHMARK(2,iit) = corr(m_max,gt_flat);
                            
                            
                            % dimensions: (mim/mic, pipeline,perfomance measure,iit)
                            PERFORMANCE(1,1:5,2,iit) = ccm.fixed;
                            PERFORMANCE(2,1:5,2,iit) = cmm.fixed;
                            PERFORMANCE(1,6,2,iit) = ccm.max;
                            PERFORMANCE(2,6,2,iit) = cmm.max;
                            PERFORMANCE(1,7,2,iit) = ccm.percent;
                            PERFORMANCE(2,7,2,iit) = cmm.percent;
                            PERFORMANCE(1,8,2,iit) = ccm.c2;
                            PERFORMANCE(2,8,2,iit) = cmm.c2;
                            
                            
                            
                            % (3) mean and std
                            for ii = 1:5
                                mean_mic.fixed(ii) = mean(mic_pca.fixed{ii}(:));
                                mean_mim.fixed(ii) = mean(mim_pca.fixed{ii}(:));
                                std_mic.fixed(ii) = std(mic_pca.fixed{ii}(:));
                                std_mim.fixed(ii) = std(mim_pca.fixed{ii}(:));
                            end
                            mean_mic.max = mean(mic_pca.max(:));
                            mean_mim.max = mean(mim_pca.max(:));
                            mean_mic.percent = mean(mic_pca.percent(:));
                            mean_mim.percent= mean(mim_pca.percent(:));
                            mean_mic.c2 = mean(mic_2(:));
                            mean_mim.c2 = mean(mim_2(:));
                            std_mic.max = std(mic_pca.max(:));
                            std_mim.max = std(mim_pca.max(:));
                            std_mic.percent = std(mic_pca.percent(:));
                            std_mim.percent= std(mim_pca.percent(:));
                            std_mic.c2 = std(mic_2(:));
                            std_mim.c2 = std(mim_2(:));
                            
                            % dimensions: (mim/mic,pipeline,perfomance measure,iit)
                            PERFORMANCE(1,1:5,3,iit) = mean_mic.fixed;
                            PERFORMANCE(2,1:5,4,iit) = std_mic.fixed;
                            PERFORMANCE(1,1:5,3,iit) = mean_mim.fixed;
                            PERFORMANCE(2,1:5,4,iit) = std_mim.fixed;
                            PERFORMANCE(1,6,[3 4],iit) = [mean_mic.max std_mic.max];
                            PERFORMANCE(2,6,[3 4],iit) = [mean_mim.max std_mim.max];
                            PERFORMANCE(1,7,[3 4],iit) = [mean_mic.percent std_mic.percent];
                            PERFORMANCE(2,7,[3 4],iit) = [mean_mim.percent std_mim.percent];
                            PERFORMANCE(1,8,[3 4],iit) = [mean_mic.c2 std_mic.c2];
                            PERFORMANCE(2,8,[3 4],iit) = [mean_mim.c2 std_mim.c2];
                            BENCHMARK([3 4],iit) = [mean(benchmark(:)) std(benchmark(:))];
                            
                            toc
                        end
                        
                        outname = sprintf('%smim_%dInts_%dRegs_0%dSNR_0%dSs.mat',DIROUT,...
                            params.iInt,params.iReg,params.isnr*10,params.iss*10);
                        save(outname,'PERFORMANCE','-v7.3')
                        
                        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
                    end
                end
            end
        end
    end
end
