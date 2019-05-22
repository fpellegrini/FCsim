function fp_megmeg_pipeline

fp_addpath

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

nit= 1000;
COH = [];
true_coh=[];

for id = 1:numel(patientID)
    
    load(sprintf('Filter_Patient%s.mat',patientID{id}));% 1D-A and true CS
    clear A
    CS = CS(1:(end-nlfp),1:(end-nlfp),:); %throw away lfp channels
    
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));    
    X = D(:,:,:);
    D_ft = ftraw(D);
    n_trials = length(D_ft.trial);
    
    id_meg_chan = 1:125;
    id_meg_chan(D.badchannels)=[];
    nmeg = numel(id_meg_chan);
    id_lfp_chan = 126:131;
    nlfp = numel(id_lfp_chan);
    
    fs = D.fsample;
    fres = 75;
    frqs = sfreqs(fres, fs);
    frqs(frqs>90) = [];
    nfreq = numel(frqs);
    
    %construct filters
    
    %leadfield
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    L1 = inverse.MEG.L;
    ns = numel(L1);    
    for is=1:ns
        L(:,is,:)= L1{is};
    end
        
    A = nan(nmeg,3,ns,nfreq);
    for ifrq = 1:nfreq      
        clear currentCS lambda CSinv
        currentCS = squeeze(CS(:,:,ifrq)); %nmeg x nmeg x nfq
        lambda = mean(diag(real(currentCS)))/100;
        CSinv=pinv(real(currentCS)+lambda * eye(size(currentCS)));
        
        for is=1:ns %iterate across nodes
            clear Lloc
            Lloc=squeeze(L(:,is,:));
            A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
        end
    end    
    
    for iit = 1:nit
        
        clear CS
        id_trials_1 = 1:n_trials;
        rng('shuffle')
        id_trials_2 = randperm(n_trials);
        CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
        CS = CS(1:(end-nlfp),1:(end-nlfp),:);
 
        %project cross spectrum to voxel space
        for ifq = 1:nfreq
            for idim = 1:3
                
                CSv = squeeze(A(:,idim,:,ifq))' * ...
                    CS(:,:,ifq) * squeeze(A(:,idim,:,ifq)); %3dim x nvox x nvox x nfreq
                
                pv = diag(squeeze(CSv));
                
                coh(idim, :, :,ifq) = squeeze(CSv) ./ sqrt(pv * pv');
                clear CSv pv
                
            end
        end
        
        %here zscore? 
        
        %aggregate voxels of one region 
        
        
        
        
    end
end


%     load(sprintf('Filter_Patient%s.mat',patientID{id}));%Filter and whole CS
%     clear CS
%     D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
%
%     X = D(:,:,:);
%     D_ft = ftraw(D);
%     n_trials = length(D_ft.trial);
%
%     fs = D.fsample;
%     fres = 75;
%     frqs = sfreqs(fres, fs);
%     frqs(frqs>90) = [];
%     nfreq = numel(frqs);
%
%     id_meg_chan = 1:125;
%     id_meg_chan(D.badchannels)=[];
%     id_trials_1 = 1:n_trials;
%
%     %get symmetric head
%     mni_pos = fp_getMNIpos(patientID{id});
%     [~, noEq] = fp_symmetric_vol(mni_pos);
%     %keep only voxels that exist in all subjects
%     A(:,noEq,:) = [];
%     A_common = A(:,voxID{id},:);
%
%     for iit = 1:nit
%         tic
%         if iit ==1
%             shuffle = 0;
%             id_trials_2 = id_trials_1;
%         else
%             shuffle =1;
%             rng('shuffle')
%             id_trials_2 = randperm(n_trials);
%         end
%
%         CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
%
%         %project cross spectrum to voxel space
%         for ifq = 1:nfreq
%             CSv(ifq,:,:) = A_common(:,:,ifq)' * CS(:,:,ifq) * A_common(:,:,ifq);
%         end
%
%         %get voxel power
%         for ifq=1:nfreq
%             pv(:,ifq) = diag(squeeze(CSv(ifq,:,:)));
%         end
%
%         %coherence
%         coh = CSv;
%         for ifreq = 1:nfreq
%             coh(ifreq, :, :) = squeeze(CSv(ifreq, :, :)) ...
%                 ./ sqrt(pv(:,ifreq)*pv(:,ifreq)');
%         end
%
%         if shuffle == 0
% %             outname = sprintf('%true_megmeg_coh_Patient%s',DIROUT, patientID{id});
% %             save(outname,'coh','-v7.3')
%             true_coh = true_coh+coh;
%         else
%             coh1(iit-1,:,:,:) = abs(imag(coh)); %nit x nfreq x nvox x nvox
%         end
%         clearvars -except patientID id nit coh1 true_coh A_common X n_trials nfreq fres n_trials id_meg_chan id_trials_1 voxID noEq
%         toc
%     end
%
%     COH = COH + coh1;
%
%     clearvars -except patientID id nit COH voxID nfreq true_coh
% end
%
% conn = fp_find_neighbours(patientID{id});
% match_conn = conn(voxID{id},voxID{id}); %is for all patients the same
% conn_s = sparse(match_conn); %nvox x nvox
% ns = size(match_conn,1);
%
% freq_conn = zeros(nfreq,nfreq);
% for ifreq = 1:nfreq-1
%     freq_conn(ifreq,ifreq+1)=1;
%     freq_conn(ifreq+1,ifreq)=1;
% end
% freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
% kron_conn = kron(conn_s,freq_conn_s); %say nvox*nfreq = nkron
% kron_conn = kron(conn_s,kron_conn);
%
% threshold = prctile(reshape(COH,1,[]),99);
%
% %true cluster
%
% clear u ind NB
% onoff_true = true_coh > threshold; %nfreq x ns
% u = onoff_true(:); %should be the same indexing like in kron_conn now; nkron x 1
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
% clu = reshape(clu,[nfreq ns ns]); %nfreq x ns x ns
%
% %save true cluster for later
% clear true_clu true_total true_sizes
% true_clu = clu;
% true_total = numel(x);
% true_sizes = x;
%
%
% %shuffled clusters
%
% onoff = COH>threshold;
% big_clusters = zeros(nit-1,nfreq,ns,ns);
%
% %find the clusters
% for iit = 1: nit-1
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
%     clu = reshape(clu,[nfreq ns ns]); %nfreq x ns
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
%         trueCoh = sum(sum(sum(true_coh(true_clu==iclus)))); %scalar
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
%
% outname = sprintf('%sp_megmeg_%s',DIROUT,abs_imag);
% save(outname,'p','threshold','true_clu','true_sizes','-v7.3')
%
%
