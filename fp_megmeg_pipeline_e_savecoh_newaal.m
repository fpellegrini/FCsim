function fp_megmeg_pipeline_e_savecoh_newaal(patientNumber,DIROUT)

fp_addpath_sabzi

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

if ~exist('DIROUT','var')
    error('Please indicate where the results should be saved.')
end

DIRLOG = '~/log/megmeg_savecoh/';
if ~exist(DIRLOG); mkdir(DIRLOG); end

ndim=2;
nit= 1000;
npcs = 5;
fres = 75;

%%
for id = 1:numel(patientID)
    fprintf('Working on subject %d. \n',id)
    logname = sprintf('%s',patientID{id});
    
    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
        tic
%%        
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
        
        %construct filters
        
        %load true CS
        load(sprintf('Filter_Patient%s_e.mat',patientID{id}));% 1D-A and true CS
        clear A
        CS = CS(1:(end-nlfp),1:(end-nlfp),:); %throw away lfp channels
        CS(:,:,1)=[];
        nfreq = size(CS,3);
        
        %leadfield
        clear L
        load(sprintf('BF_Patient%s.mat',patientID{id}));
        L = fp_get_lf(inverse);
        ns_org = size(L,2);
        
        %get rois
        clear mni_pos label code roi_id u_roi_id csroi
        mni_pos = fp_getMNIpos(patientID{id});
        for ii = 1: ns_org
            [~,~,roi_id(ii)]=fp_get_mni_anatomy_new(mni_pos(ii,:));
        end
        u_roi_id = sort(unique(roi_id));
        nroi = numel(u_roi_id);
        
        %get rid of white voxels
        L(:,roi_id==0,:)=[];
        ns = size(L,2);
        
        %construct beamformer
        A = squeeze(mkfilt_eloreta_v2(L));
        A = permute(A,[1, 3, 2]);
        
        % true coherence
        coh_roi = nan(nroi-1,nroi-1,npcs,npcs,nfreq);
        
        clear cA A_ CSn v5 cseig
        cA = A;
        A_ = reshape(cA, [nmeg, ndim*ns]);
        for ifq = 1:nfreq
            clear pv
            CSv(ifq,:,:) = A_' * CS(:,:,ifq) * A_;
%             pv = fp_project_power(CS(:,:,ifq),A_);
        end       
        
        %zscoring
        ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 1))))./diag(squeeze(sum(real(CSv), 1)))));
        for ifreq = 1:nfreq
            CSz(ifreq,:, :) = ZS'*squeeze(CSv(ifreq,:, :))*ZS;
        end
        
        %region pca
        CSs = squeeze(sum(CSz,1));       
        is = 1;
        for iroi = 2:nroi
            clear v cCS cns iid
            
            iid = is: is+ (sum(roi_id == u_roi_id(iroi)))*ndim-1;
            cCS = CSs(iid,iid);
            [v, ~, ~] = eig(real(cCS));
            
            v5{iroi} = v(:,1:npcs); %nregionvoxels*2 x npcs
            is = is+length(iid);
        end
        clear CSs
        
        %apply the pca-filters to the cs
        for ifq = 1:nfreq
            
            clear cseig
            kr=1;
            for kroi = 2:nroi
                clear kid
                
                kid = kr: kr+ (sum(roi_id == u_roi_id(kroi)))*ndim-1;
                jr=1;
                for jroi =2: nroi
                    clear cCS jid
                    jid = jr:jr+ (sum(roi_id == u_roi_id(jroi)))*ndim-1;
                    cCS = squeeze(CSz(ifq,kid,jid));
                    
                    cseig(kroi-1,jroi-1,:,:) = v5{kroi}' * cCS * v5{jroi};
                    
                    jr=jr+length(jid);
                end
                
                kr=kr+length(kid);
            end
            
            %divide by power
            for ifc=1:npcs
                for jfc =1:npcs
                    clear proi
                    proi = squeeze(real(diag(cseig(:,:,ifc,jfc))));
                    coh_roi(:,:,ifc,jfc,ifq)= cseig(:,:,ifc,jfc)./sqrt(proi' * proi);
                end
            end
            
        end
        
        clear true_coh
        for ii=1:nroi-1
            for jj= 1: nroi-1
                for ifq = 1: nfreq
                    true_coh(ii,jj,ifq) =  sum(sum(squeeze(abs(imag(coh_roi(ii,jj,:,:,ifq))))));
                end
            end
        end
        
 %%       
        
        % calculate coherence for permutations
        
        for iit = 1:nit
            fprintf('Working on iteration %d. \n',iit)
            
            %cross spectrum
            
            clear CS coh
            id_trials_1 = 1:n_trials;
            rng('shuffle')
            id_trials_2 = randperm(n_trials);
            CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
            
            clear A_ coh_roi
            coh_roi = nan(nroi-1,nroi-1,npcs,npcs,nfreq);
            A_ = reshape(A, [nmeg, ndim*ns]);
            clear CSv pv CSz cseig
            
            for ifq = 1:nfreq              
                CSv(ifq,:,:) = A_' * CS(:,:,ifq) * A_;
%                 pv = fp_project_power(CS(:,:,ifq),A_);
%                 CSz = CSv ./ sqrt(pv * pv');
            end
            
            %zscoring
            ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 1))))./diag(squeeze(sum(real(CSv), 1)))));
            
            for ifreq = 1:nfreq
                
                CSz = ZS'*squeeze(CSv(ifreq,:, :))*ZS;
                               
                %apply the filters to the cs
                kr=1;
                for kroi = 2:nroi
                    clear kid
                    
                    kid = kr: kr+ (sum(roi_id == u_roi_id(kroi)))*ndim-1;
                    jr=1;
                    for jroi =2: nroi
                        clear cCS jid
                        jid = jr:jr+ (sum(roi_id == u_roi_id(jroi)))*ndim-1;
                        cCS = CSz(kid,jid);
                        
                        cseig(kroi-1,jroi-1,:,:) = v5{kroi}' * cCS * v5{jroi};
                        
                        jr=jr+length(jid);
                    end
                    
                    kr=kr+length(kid);
                end
                
                %divide by power
                for ifc=1:npcs
                    for jfc =1:npcs
                        clear proi
                        proi = squeeze(real(diag(cseig(:,:,ifc,jfc))));
                        coh_roi(:,:,ifc,jfc,ifq)= cseig(:,:,ifc,jfc)./sqrt(proi' * proi);
                    end
                end
            end
            
            clear coh
            
            for ii=1:nroi-1
                for jj= 1: nroi-1
                    for ifq = 1: nfreq
                        coh(ii,jj,ifq) =  sum(sum(squeeze(abs(imag(coh_roi(ii,jj,:,:,ifq))))));
                    end
                end
            end
            
            COH_sub(iit,:,:,:) = coh;
        end
        
        outname = sprintf('%sroi_coh_sub%s',DIROUT,patientID{id});
        save(outname,'COH_sub','true_coh','-v7.3')
        
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
        clearvars -except id patientID nit npcs nfreq ndim nroi DIRLOG DIROUT
        toc
    end
end

%
% %%
% %neighbourhood
% load('roi_conn.mat')
% roiconn_s = sparse(roi_conn);
% freq_conn = fp_get_freq_conn(nfreq);
% freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
% kron_conn = kron(roiconn_s,freq_conn_s); %nroi*nfreq = nkron
% kron_conn = kron(roiconn_s,kron_conn);
%
% threshold = prctile(COH(:),99.9);
%
% %%
% %true cluster
% clear onoff
% onoff = TRUE_COH>threshold;
%
% %find the clusters
% clear u ind NB
% u = onoff(:); %should be the same indexing like in kron_conn now; nkron x 1
%
% ind = find(u==1); %remember indeces of super-threshold coherences
% NB = kron_conn;
% NB(u==0,:)=[]; %pass the neighbourhood structure only for the super-threshold voxels
% NB(:,u==0)=[];
%
% %components assigns every voxel to a cluster, even if this means that every voxel is its own cluster
% clear ci x clu
% [ci, x] = components(NB); %x is the histogram of the clusters
% clu = zeros(size(kron_conn,1),1);%refill with sub-threshold voxels
% clu(ind)= ci; %nkron x 1
% clu = reshape(clu,[nfreq nroi-1 nroi-1]);
%
% %save true cluster for later
% clear true_clu_pos true_total_pos
% true_clu = clu;
% true_total = numel(x);
%
%
%
% %% shuffled clusters
%
% clear onoff
% onoff = COH>threshold;
% big_clusters = zeros(nit,nfreq,nroi-1,nroi-1);
%
% %find the clusters
% for iit = 1: nit
%
%     clear onoff_temp u ind NB
%     onoff_temp = squeeze(onoff(iit,:,:)); %nfreq x ns
%     u = onoff_temp(:); %should be the same indexing like in kron_conn now; nkron x 1
%
%     ind = find(u==1); %remember indeces of super-threshold coherences
%     NB = kron_conn;
%     NB(u==0,:)=[]; %pass the neighbourhood structure only for the super-threshold voxels
%     NB(:,u==0)=[];
%
%     %components assigns every voxel to a cluster, even if this means that every voxel is its own cluster
%     clear ci x clu
%     [ci, x] = components(NB); %x is the histogram of the clusters
%     clu = zeros(size(kron_conn,1),1);%refill with sub-threshold voxels
%     clu(ind)= ci; %nkron x 1
%     clu = reshape(clu,[nfreq nroi-1 nroi-1]);
%
%     if numel(x)>0
%         big_clu_id = find(x==max(x));
%         big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
%         big_clusters(iit,:,:,:) = (clu == big_clu_id);
%     end
%
% end
%
% %compare not only cluster size but also magnitude of coherence within
% %the relevant cluster
%
% clear a
% a = zeros(size(COH)); %only shuffled clusters
% a(big_clusters==1) = COH(big_clusters==1);
% shufCoh = squeeze(sum(sum(sum(a,2),3),4));
%
% if true_total>0 %when at least one true cluster exists
%     for iclus = 1:true_total
%         clear trueCoh temp
%         trueCoh = sum(sum(sum(TRUE_COH(true_clu==iclus)))); %scalar
%         p(iclus) = sum(shufCoh>trueCoh)/numel(shufCoh);
%     end
%
% elseif sum(shufCoh)== 0  %when no cluster was found it any iteration
%     p= nan;
%
% else %when only in shuffled conditions clusters were found
%     clear trueCoh
%     trueCoh = 0;
%     p = sum(shufCoh>trueCoh)/numel(shufCoh);
% end
%
% %%
% outname = sprintf('%sp_megmeg_%s',DIROUT,abs_imag);
% save(outname,'p','threshold','true_clu','-v7.3')
%

