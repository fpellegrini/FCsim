function fp_megmeg_pipeline_mim(patientNumber,DIROUT)

fp_addpath_sabzi

if ~exist(DIROUT); mkdir(DIROUT); end

DIRLOG = '~/log/mim90/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

if ~exist('DIROUT','var')
    error('Please indicate where the results should be saved.')
end

nit= 1000;

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
        id_trials = 1:n_trials;
        
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
        
        %CS
        CS = fp_tsdata_to_cpsd(X,fres,'WELCH',[id_meg_chan], [id_meg_chan], id_trials, id_trials);
        CS(:,:,[1 47:end])=[];
        nfreq = size(CS,3);
        
        %leadfield
        clear L
        load(sprintf('BF_Patient%s.mat',patientID{id}));
        L = fp_get_lf(inverse);
        [~, ns_org, ni] = size(L);
        
        %get rois
        clear mni_pos label code roi_id u_roi_id csroi
        mni_pos = fp_getMNIpos(patientID{id});
        for ii = 1: ns_org
            [label{ii},code{ii},roi_id(ii),partner_rois]=fp_get_mni_anatomy_new(mni_pos(ii,:));
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
        
        % true coherence
        clear P V A2 
        for aroi = 1:nroi
            
            %project to source level
            clear A_ CSv
            A_ = A(:, :,roi_id == aroi,:);
            nvoxroi = size(A_,3);
            A2{aroi} = reshape(A_, [nmeg, ni*nvoxroi, nfreq]);
            
            for ifq = 1:nfreq
                CSv(:,:,ifq) = squeeze(A2{aroi}(:,:,ifq))' * CS(:,:,ifq)...
                    * squeeze(A2{aroi}(:,:,ifq));
            end
            
            %zscoring
            clear CSz
            ZS{aroi} = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
            for ifreq = 1:nfreq
                CSz(ifreq,:, :) = ZS{aroi}'*squeeze(CSv(:,:, ifreq))*ZS{aroi};
            end
            
            %region pca
            clear CSs v v5 in V_ D_
            CSs = squeeze(sum(CSz,1)); %covariance           
            [V_, D_] = eig(real(CSs));
            [D_, in] = sort(real(diag(D_)), 'descend');
            % variance explained
            vx_ = cumsum(D_)./sum(D_);
            invx = 1:min(length(vx_), nmeg);
            npcs(aroi) = min(find(vx_>0.9));
            
            V{aroi} = V_(:,in);

        end
        
        %makes sure that rois have the same npcs in both hemispheres 
        croi = 1;
        for aroi = 1:nroi 
            if ~isnan(partner_rois(2,aroi))
                npcs(aroi) = max(npcs(aroi),npcs(partner_rois(2,aroi)));
            end
            V{aroi} = V{aroi}(:, 1:npcs(aroi));
            for ifq = 1:nfreq
                P(:, croi:croi+npcs(aroi)-1,ifq) = A2{aroi}(:,:,ifq) * ZS{aroi} * real(V{aroi});
            end
            croi = croi +npcs(aroi);
        end
        clear ZS
                
        %apply all filters
        CSroi = [];
        for ifreq = 1:nfreq
            CSroi(:, :, ifreq) = reshape(P(:,:,ifreq), nmeg, [])'*CS(:, :, ifreq)...
                *reshape(P(:,:,ifreq), nmeg, []);
        end
        
        %divide by power to obtain coherence
        clear Cohroi
        for ifreq = 1: nfreq
            clear pow
            pow = real(diag(CSroi(:,:,ifreq)));
            Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
        end
        
        
        %
        %mim and mic
        clear mic mim
        [mic,mim] =  fp_mim(Cohroi,npcs);
        
        MIC_TRUE(id,:,:,:) = mic;
        MIM_TRUE(id,:,:,:) = mim;
    %%    
        % permutations
        
        for iit = 1:nit
            fprintf('Working on iteration %d. \n',iit)
            
            tic
            %cross spectrum
            clear CS coh
            id_trials_1 = 1:n_trials;
            id_trials_2 = randperm(n_trials);
            CS = fp_tsdata_to_cpsd(X,fres,'WELCH',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
            
            clear P_shuf
            croi = 1; 
            for aroi = 1:nroi
                
               %project to source level
               for ifq = 1:nfreq
                    CSv(:,:,ifq) = squeeze(A2{aroi}(:,:,ifq))' * CS(:,:,ifq)...
                        * squeeze(A2{aroi}(:,:,ifq));
                end
                
                %zscoring
                clear CSz ZS
                ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))...
                    ./diag(squeeze(sum(real(CSv), 3)))));
                for ifreq = 1:nfreq
                    CSz(ifreq,:, :) = ZS'*squeeze(CSv(:,:, ifreq))*ZS;
                end
                
                for ifq = 1:nfreq
                    P_shuf(:, croi:croi+npcs(aroi)-1,ifq) = A2{aroi}(:,:,ifq) * ZS * real(V{aroi});
                end
                croi = croi +npcs(aroi);
            end
            
            %apply all filters
            CSroi = [];
            for ifreq = 1:nfreq
                CSroi(:, :, ifreq) = reshape(P_shuf(:,:,ifreq), nmeg, [])'*CS(:, :, ifreq)...
                    *reshape(P_shuf(:,:,ifreq), nmeg, []);
            end
            
            %divide by power to obtain coherence
            clear Cohroi
            for ifreq = 1: nfreq
                clear pow
                pow = real(diag(CSroi(:,:,ifreq)));
                Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
            end
            
            %mim and mic
            clear mic mim
            [mic,mim] =  fp_mim(Cohroi,npcs);
            
            MIC_SHUF(iit,id,:,:,:) = mic;
            MIM_SHUF(iit,id,:,:,:) = mim;
            toc
        end
        
        outname = sprintf('%sroi_MIM90_sub%s',DIROUT,patientID{id});
        save(outname,'MIC_TRUE','MIM_TRUE','MIC_SHUF','MIM_SHUF','-v7.3')
        
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
    end
end
