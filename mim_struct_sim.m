%todo:
%1) evt replace tsdata with Guido's code?

%%

DIROUT = './';
DIRLOG = './log/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

rng('shuffle')
varyParam = 1%:8;
nit = 1%500;
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
    
    
    for iInt = nInteractions
        for iReg = nRegionInts
            for isnr = SNR
                for iss = noise_mix
                    for ilag = nlag
                        for ifilt = filtertype
                            for ihemi = hemisym
                                
                                %create logfile for parallization
                                logname = sprintf('iInt%d_iReg%d_snr0%d_iss0%d_lag%d_filt%s_hemisym%d'...
                                    ,iInt,iReg,isnr*10,iss*10, ilag,ifilt,ihemi);
                                
                                if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
                                    eval(sprintf('!touch %s%s_work',DIRLOG,logname))
                                    fprintf('Working on %s. \n',logname)
                                    
                                    clear params PERFORMANCE
                                    params.iInt = iInt;
                                    params.iReg = iReg;
                                    params.isnr = isnr;
                                    params.iss = iss;
                                    params.ilag = ilag;
                                    params.ifilt = ifilt;
                                    params.ihemi = ihemi;
                                    
                                    % dimensions: (mim/mic, pipeline,perfomance measure,iit)
                                    PERFORMANCE = zeros(2,8,2,nit);
                                    
                                    %in each iteration, a new signal with new
                                    %interacting sources and voxels is generated
                                    for iit = 1: nit
                                        
                                        fprintf('Working on iteration %d. \n',iit)
                                        %% signal generation
                                        
                                        clearvars -except DIROUT DIRLOG mm_gt mc_gt bmm_gt bmc_gt GT MIM MIC ...
                                            params iit nit nInteractions nRegionInts SNR noise_mix ...
                                            nlag filtertype hemisym ip varyParam...
                                            fres n_trials PERFORMANCE iss iInt iReg isnr ilag ifilt ihemi
                                        
                                        
                                        % ROI labels
                                        % In D.sub_ind_roi, there are the randomly
                                        % selected voxels of each region
                                        fprintf('Getting atlas positions... \n')
                                        tic
                                        D = fp_get_Desikan(params.iReg);
                                        toc
                                        
                                        fprintf('Signal generation... \n')
                                        tic
                                        %signal generation
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
                                        gt_save = gt;
                                        mic_save = mic;
                                        mim_save = mim;
                                        
                                        gt.mic = sum(gt.mic,3);
                                        gt.mim = sum(gt.mim,3);
                                        
                                        for ii =1:5
                                            mic.fixed{ii} = sum(mic.fixed{ii},3);
                                            mim.fixed{ii} = sum(mim.fixed{ii},3);
                                        end
                                        mic.max = sum(mic.max,3);
                                        mim.max = sum(mim.max,3);
                                        mic.percent = sum(mic.percent,3);
                                        mim.percent = sum(mim.percent,3);
                                        mic.case2 = sum(mic.case2,3);
                                        mim.case2 = sum(mim.case2,3);
                                        mic.baseline = sum(mic.baseline,3);
                                        mim.baseline = sum(mim.baseline,3);
                                        
                                        
                                        % (1)correlation mim/mic and ground truth
                                        for ii = 1:5
                                            PERFORMANCE(1,ii,1,iit) = corr(mic.fixed{ii}(:),gt.mic(:));
                                            PERFORMANCE(2,ii,1,iit) = corr(mim.fixed{ii}(:),gt.mim(:));
                                        end
                                        PERFORMANCE(1,6,1,iit) = corr(mic.max(:),gt.mic(:));
                                        PERFORMANCE(1,7,1,iit) = corr(mic.percent(:),gt.mic(:));
                                        PERFORMANCE(1,8,1,iit) = corr(mic.case2(:),gt.mic(:));
                                        PERFORMANCE(2,6,1,iit) = corr(mim.max(:),gt.mim(:));
                                        PERFORMANCE(2,7,1,iit) = corr(mim.percent(:),gt.mim(:));
                                        PERFORMANCE(2,8,1,iit) = corr(mim.case2(:),gt.mim(:));
                                        
                                        %baseline
                                        BASELINE(1,1,iit) = corr(mic.baseline(:),gt.mic(:));
                                        BASELINE(2,1,iit) = corr(mim.baseline(:),gt.mim(:));
                                        
                                        
                                        % (2) correlation maxima of mim/mic and ground truth
                                        for ii = 1:5
                                            %mic
                                            clear m_max
                                            m_max = fp_get_nmaxima(mic.fixed{ii},params.iInt*2);
                                            PERFORMANCE(1,ii,2,iit) = corr(m_max(:),gt.mic(:));
                                            
                                            %mim
                                            clear m_max
                                            m_max = fp_get_nmaxima(mim.fixed{ii},params.iInt*2);
                                            PERFORMANCE(2,ii,2,iit) = corr(m_max(:),gt.mim(:));
                                        end
                                        
                                        %max pipeline
                                        %mic
                                        clear m_max
                                        m_max = fp_get_nmaxima(mic.max,params.iInt*2);
                                        PERFORMANCE(1,6,2,iit) = corr(m_max(:),gt.mic(:));
                                        %mim
                                        clear m_max
                                        m_max = fp_get_nmaxima(mim.max,params.iInt*2);
                                        PERFORMANCE(2,6,2,iit) = corr(m_max(:),gt.mim(:));
                                        
                                        %percent pipeline
                                        %mic
                                        clear m_max
                                        m_max = fp_get_nmaxima(mic.percent,params.iInt*2);
                                        PERFORMANCE(1,7,2,iit) = corr(m_max(:),gt.mic(:));
                                        %mim
                                        clear m_max
                                        m_max = fp_get_nmaxima(mim.percent,params.iInt*2);
                                        PERFORMANCE(2,7,2,iit) = corr(m_max(:),gt.mim(:));
                                        
                                        %case2 pipeline
                                        %mic
                                        clear m_max
                                        m_max = fp_get_nmaxima(mic.case2,params.iInt*2);
                                        PERFORMANCE(1,8,2,iit) = corr(m_max(:),gt.mic(:));
                                        %mim
                                        clear m_max
                                        m_max = fp_get_nmaxima(mim.case2,params.iInt*2);
                                        PERFORMANCE(2,8,2,iit) = corr(m_max(:),gt.mim(:));
                                        
                                        %baseline
                                        %mic
                                        clear m_max
                                        m_max = fp_get_nmaxima(mic.baseline,params.iInt*2);
                                        BASELINE(1,2,iit) = corr(m_max(:),gt.mic(:));
                                        %mim
                                        clear m_max
                                        m_max = fp_get_nmaxima(mim.baseline,params.iInt*2);
                                        BASELINE(2,2,iit) = corr(m_max(:),gt.mim(:));
                                        toc
                                        
                                    end
                                    
                                    fprintf('Saving... \n')
                                    tic
                                    outname = sprintf('%smim_%s.mat',DIROUT,logname);
                                    save(outname,'PERFORMANCE','BASELINE','-v7.3')
                                    toc
                                    
                                    eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
