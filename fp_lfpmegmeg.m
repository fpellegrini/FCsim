function fp_lfpmegmeg

fp_addpath_sabzi

DIRLOG = '~/log/lfpmegmeg1/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

DIROUT = '~/data/lfpmegmeg1/';
if ~exist(DIROUT); mkdir(DIROUT); 
end

%%
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
nit= 100;

%this is Precentral left and right, SMA left and right, parietal left and
%right, cerebellum, and pallidum
myrois = [1 15 3 17 8 22 13 14]; 

fixed_npcs = 5;


%%
for id = 1:numel(patientID)
    clearvars -except DIRLOG DIROUT patientID nit myrois id fixed_npcs
    
    logname = sprintf('%s',patientID{id});
    
    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
        fprintf('Working on subject %d. \n',id)
        
        %load data
        clear X V ZS A2
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
        z = exp(-i*pi*(frqs./(fs/2)))';
        
        %construct filters
        
        %CS
        CS = fp_tsdata_to_cpsd(X,fres,'WELCH',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan], id_trials, id_trials);
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
            [~,~,roi_id(ii),partner_rois]=fp_get_mni_anatomy_new(mni_pos(ii,:));
        end
        %get only my rois 
        roi_id_save = roi_id; 
        roi_id(~ismember(roi_id, myrois))=0; 
        u_roi_id = sort(unique(roi_id));
        nroi = numel(u_roi_id)-1;       
        %get rid of all other rois 
        L(:,roi_id ==0, :) = [];
        roi_id(roi_id==0) = []; 
        partner_rois = partner_rois(:,myrois);
        nvox = size(L,2);
        
        %construct lcmv filter (lcmv because we also want to calculate the
        %GC later) 
        cCS = sum(CS(1:nmeg,1:nmeg,:),3);
        reg = 0.05*trace(cCS)/length(cCS);
        Cr = cCS + reg*eye(size(cCS,1));        
        [~, A] = lcmv_meg(Cr, L, struct('alpha', 0, 'onedim', 0));
        A = permute(A,[1, 3, 2]);
        
        % true coherence
        clear P V A2 
        for aroi = 1:nroi
            
            %project to source level
            clear A_ CSv
            A_ = A(:, :,roi_id == myrois(aroi));
            nvoxroi = size(A_,3);
            A2{aroi} = reshape(A_, [nmeg, ni*nvoxroi]);
            
            for ifq = 1:nfreq
                CSv(:,:,ifq) = squeeze(A2{aroi})' * CS(1:nmeg,1:nmeg,ifq)...
                    * squeeze(A2{aroi});
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
%             vx_ = cumsum(D_)./sum(D_);
%             invx = 1:min(length(vx_), nmeg);
            npcs(aroi) = fixed_npcs; %min(find(vx_>0.9));
            
            V{aroi} = V_(:,in);

        end
        
        %makes sure that rois have the same npcs in both hemispheres 
        croi = 1;
        for aroi = 1:nroi 
%             if ~isnan(partner_rois(2,aroi))
%                 npcs(aroi) = max(npcs(aroi),npcs(partner_rois(2,aroi)));
%             end
            V{aroi} = V{aroi}(:, 1:npcs(aroi));
            P(:, croi:croi+npcs(aroi)-1) = A2{aroi} * ZS{aroi} * real(V{aroi});
            croi = croi +npcs(aroi);
        end
        clear ZS
                
        %apply all filters - meg and lfp
        croi = croi-1;
        CSroi = zeros(croi+nlfp,croi+nlfp,nfreq);
        for ifreq = 1:nfreq
            CSroi(1:croi, 1:croi, ifreq) = P'*CS(1:nmeg, 1:nmeg, ifreq)*P; %megmeg
            CSroi(1:croi,croi+1:end,ifreq) = P'* CS(1:nmeg,nmeg+1:end,ifreq); %meglfp
            CSroi(croi+1:end,1:croi,ifreq) = CS(nmeg+1:end,1:nmeg,ifreq) * P; %lfpmeg
            CSroi(croi+1:end,croi+1:end,ifreq) = CS(nmeg+1:end,nmeg+1:end,ifreq); %lfplfp
        end
        
        %divide by power to obtain coherence
        clear Cohroi
        for ifreq = 1: nfreq
            clear pow
            pow = real(diag(CSroi(:,:,ifreq)));
            Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
        end
       
        %npcs for lfp channels are 3 (3 lfps on right and 3 lfps on left
        %side)
        npcs = [npcs nlfp/2 nlfp/2];
        
        %mim and mic
        clear mic mim
        [mic,mim] =  fp_mim(Cohroi,npcs);
        
        MIC_TRUE = mic;
        MIM_TRUE = mim;
        
        %granger causality
        [GC_TRUE, TRGC_TRUE, DIFFGC_TRUE] = fp_gc(CSroi,npcs,z); 
        
        
        
        %% permutations
        
        for iit = 1:nit
            fprintf('Working on iteration %d. \n',iit)
            
            tic
            %cross spectrum
            clear CS coh
            id_trials_1 = 1:n_trials;
            id_trials_2 = randperm(n_trials);
            CS = fp_tsdata_to_cpsd(X,fres,'WELCH',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan], id_trials_1, id_trials_2);
            CS(:,:,[1 47:end])=[];
             
            clear P_shuf
            croi = 1; 
            for aroi = 1:nroi
                
               %project to source level
               clear CSv
               for ifq = 1:nfreq
                    CSv(:,:,ifq) = squeeze(A2{aroi})' * CS(1:nmeg,1:nmeg,ifq)...
                        * squeeze(A2{aroi});
                end
                
                %zscoring
                clear CSz ZS a
                a = diag(squeeze(sum(real(CSv), 3)));    
                ZS = diag(sqrt(mean(a)./a));                    
                P_shuf(:, croi:croi+npcs(aroi)-1) = A2{aroi} *ZS* real(V{aroi});
                croi = croi +npcs(aroi);
            end

            %apply all filters - meg and lfp
            croi = croi-1;
            CSroi = zeros(croi+nlfp,croi+nlfp,nfreq);
            for ifreq = 1:nfreq
                CSroi(1:croi, 1:croi, ifreq) = P_shuf'*CS(1:nmeg, 1:nmeg, ifreq)* P_shuf; %megmeg
                CSroi(1:croi,croi+1:end,ifreq) = P_shuf'* CS(1:nmeg,nmeg+1:end,ifreq); %meglfp
                CSroi(croi+1:end,1:croi,ifreq) = CS(nmeg+1:end,1:nmeg,ifreq) * P_shuf; %lfpmeg
                CSroi(croi+1:end,croi+1:end,ifreq) = CS(nmeg+1:end,nmeg+1:end,ifreq); %lfplfp
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
            
            MIC_SHUF(iit,:,:,:) = mic;
            MIM_SHUF(iit,:,:,:) = mim;
            toc
        end
        
        outname = sprintf('%sMIM_GC_sub%s',DIROUT,patientID{id});
        save(outname,'MIC_TRUE','MIM_TRUE','MIC_SHUF','MIM_SHUF','GC_TRUE',...
            'TRGC_TRUE','DIFFGC_TRUE','-v7.3')
        
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
    end
end
