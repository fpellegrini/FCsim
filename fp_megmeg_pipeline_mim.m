function fp_megmeg_pipeline_mim(patientNumber,DIROUT,DIRLOG)

fp_addpath_sabzi

if ~exist(DIROUT); mkdir(DIROUT); end

if ~exist(DIRLOG); mkdir(DIRLOG); end

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

if ~exist('DIROUT','var')
    error('Please indicate where the results should be saved.')
end

nit= 2;
npcs = 2;
COH = zeros(nit,117,117,46);
TRUE_COH = zeros(117,117,46);
regu=.000001;

thetaband = 2:4; %same as 4 to 8 Hz
alphaband = 4:6; % 8 to 12 Hz
betaband = 6:15; % 12 to 30 Hz
gammaband = 16:45; %32 to 90 Hz;

%%
for id = 1:numel(patientID)
    logname = sprintf('%s',patientID{id});
    
    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
        fprintf('Working on subject %d. \n',id)
        
        %load data
        clear X
        D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
        X = D(:,:,:);
        D_ft = ftraw(D);
        n_trials = length(D_ft.trial);
        
        %channel IDs
        clear id_meg_chan id_lfp_chan
        id_meg_chan = 1:125;
        id_meg_chan(D.badchannels)=[];
        nmeg = numel(id_meg_chan);
        id_lfp_chan = 126:131;
        nlfp = numel(id_lfp_chan);
        
        %scaling
        X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./(10^6);
        
        %frequency parameters
        fs = D.fsample;
        fres = 75;
        frqs = sfreqs(fres, fs);
        frqs(frqs>90) = [];
        nfreq = numel(frqs);
        
        %construct filters
        
        %load true CS
        load(sprintf('Filter_Patient%s_e.mat',patientID{id}));% 1D-A and true CS
        clear A
        CS = CS(1:(end-nlfp),1:(end-nlfp),:); %throw away lfp channels
        
        %leadfield
        clear L
        load(sprintf('BF_Patient%s.mat',patientID{id}));
        L = fp_get_lf(inverse);
        [~, ns_org, ni] = size(L);
        
        %get rois
        clear mni_pos label code roi_id u_roi_id csroi
        mni_pos = fp_getMNIpos(patientID{id});
        for ii = 1: ns_org
            [label{ii},code{ii},roi_id(ii)]=fp_get_mni_anatomy_new(mni_pos(ii,:));
        end
        u_roi_id = sort(unique(roi_id));
        nroi = numel(u_roi_id)-1;
        
        %get rid of white voxels
        %     L(:,roi_id==0,:)=[];
        nvox = size(L,2);
        
        %construct beamformer
        A = nan(nmeg,ni,nvox,nfreq);
        for ifrq = 1:nfreq
            clear cCS lambda CSinv
            cCS = squeeze(CS(:,:,ifrq));
            lambda = mean(diag(real(cCS)))/100;
            
            CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
            
            for is=1:nvox %iterate across nodes
                clear Lloc
                Lloc=squeeze(L(:,is,:));
                A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
            end
        end
        
        %% true coherence
        clear P
        
        for aroi = 1:nroi
            
            %project to source level
            clear A_ CSv A2
            A_ = A(:, :,roi_id == aroi,:);
            nvoxroi = size(A_,3);
            A2 = reshape(A_, [nmeg, ni*nvoxroi, nfreq]);
            
            for ifq = 1:nfreq
                CSv(:,:,ifq) = squeeze(A2(:,:,ifq))' * CS(:,:,ifq)...
                    * squeeze(A2(:,:,ifq));
            end
            
            %zscoring
            clear ZS CSz
            ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
            for ifreq = 1:nfreq
                CSz(ifreq,:, :) = ZS'*squeeze(CSv(:,:, ifreq))*ZS;
            end
            
            %region pca
            clear CSs v v5 in V_ D_
            CSs = squeeze(sum(CSz,1)); %covariance
            
            [V_, D_] = eig(real(CSs));
            [D_, in] = sort(real(diag(D_)), 'descend');
            
            V{aroi} = V_(:,in(1:npcs)); %nregionvoxels*2 x npcs
            
            %     %concatenate filters
            for ifq = 1:nfreq
                P(:, :, aroi,ifq) = A2(:,:,ifq) * ZS * real(V{aroi});
            end
        end
        
        %apply all filters
        CSroi = [];
        for ifreq = 1:nfreq
            CSroi(:, :, ifreq) = reshape(P(:,:,:,ifreq), nmeg, [])'*CS(:, :, ifreq)...
                *reshape(P(:,:,:,ifreq), nmeg, []);
        end
        
        %divide by power to obtain coherence
        clear Cohroi
        for ifreq = 1: nfreq
            clear pow
            pow = real(diag(CSroi(:,:,ifreq)));
            Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
        end
        
        %mim and mic
        clear mim1 mic1 mic mim
        [mic1,mim1] =  fp_mim(Cohroi);
        
        mic(:,:,1) = mean(mic1(:,:,thetaband),3);
        mic(:,:,2) = mean(mic1(:,:,alphaband),3);
        mic(:,:,3) = mean(mic1(:,:,betaband),3);
        mic(:,:,4) = mean(mic1(:,:,gammaband),3);
        mim(:,:,1) = mean(mim1(:,:,thetaband),3);
        mim(:,:,2) = mean(mim1(:,:,alphaband),3);
        mim(:,:,3) = mean(mim1(:,:,betaband),3);
        mim(:,:,4) = mean(mim1(:,:,gammaband),3);
        
        MIC_TRUE(id,:,:,:) = mic;
        MIM_TRUE(id,:,:,:) = mim;
        
        %% permutations
        
        for iit = 1:nit
            
            %cross spectrum
            clear CS coh
            id_trials_1 = 1:n_trials;
            rng('shuffle')
            id_trials_2 = randperm(n_trials);
            CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
            
            clear P_shuf
            
            for aroi = 1:nroi
                
                %project to source level
                clear A_ CSv A2
                A_ = A(:, :,roi_id == aroi,:);
                nvoxroi = size(A_,3);
                A2 = reshape(A_, [nmeg, ni*nvoxroi, nfreq]);
                
                for ifq = 1:nfreq
                    CSv(:,:,ifq) = squeeze(A2(:,:,ifq))' * CS(:,:,ifq)...
                        * squeeze(A2(:,:,ifq));
                end
                
                %zscoring
                clear ZS CSz
                ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
                for ifreq = 1:nfreq
                    CSz(ifreq,:, :) = ZS'*squeeze(CSv(:,:, ifreq))*ZS;
                end
                
                
                %     %concatenate filters
                for ifq = 1:nfreq
                    P_shuf(:, :, aroi,ifq) = A2(:,:,ifq) * ZS * real(V{aroi});
                end
            end
            
            %apply all filters
            CSroi = [];
            for ifreq = 1:nfreq
                CSroi(:, :, ifreq) = reshape(P_shuf(:,:,:,ifreq), nmeg, [])'*CS(:, :, ifreq)...
                    *reshape(P_shuf(:,:,:,ifreq), nmeg, []);
            end
            
            %divide by power to obtain coherence
            clear Cohroi
            for ifreq = 1: nfreq
                clear pow
                pow = real(diag(CSroi(:,:,ifreq)));
                Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
            end
            
            %mim and mic
            clear mim1 mic1 mic mim
            [mic1,mim1] =  fp_mim(Cohroi);
            
            mic(:,:,1) = mean(mic1(:,:,thetaband),3);
            mic(:,:,2) = mean(mic1(:,:,alphaband),3);
            mic(:,:,3) = mean(mic1(:,:,betaband),3);
            mic(:,:,4) = mean(mic1(:,:,gammaband),3);
            mim(:,:,1) = mean(mim1(:,:,thetaband),3);
            mim(:,:,2) = mean(mim1(:,:,alphaband),3);
            mim(:,:,3) = mean(mim1(:,:,betaband),3);
            mim(:,:,4) = mean(mim1(:,:,gammaband),3);
            
            MIC_SHUF(iit,id,:,:,:) = mic;
            MIM_SHUF(iit,id,:,:,:) = mim;
        end
        
        outname = sprintf('%sroi_MIM_sub%s',DIROUT,patientID{id});
        save(outname,'MIC_TRUE','MIM_TRUE','MIC_SHUF','MIM_SHUF','-v7.3')
        
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
    end
end
