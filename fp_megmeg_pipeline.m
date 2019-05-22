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
    
    %calculate coherence 
    
    for iit = 1:nit
        
        %cross spectrum
        clear CS
        id_trials_1 = 1:n_trials;
        rng('shuffle')
        id_trials_2 = randperm(n_trials);
        CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
        CS = CS(1:(end-nlfp),1:(end-nlfp),:);
 
        %project cross spectrum to voxel space and get power and coherence 
        for ifq = 1:nfreq
            for idim = 1:3
                
                CSv = squeeze(A(:,idim,:,ifq))' * ...
                    CS(:,:,ifq) * squeeze(A(:,idim,:,ifq)); %3dim x nvox x nvox x nfreq
                
                pv = diag(squeeze(CSv));
                
                coh(idim, :, :,ifq) = squeeze(CSv) ./ sqrt(pv * pv');
                clear CSv pv
                
            end
        end
        
        %%%%%%here zscore? 
        
        %aggregate voxels of one region 
        
        clear mni_pos label code roi_id u_roi_id
        mni_pos = fp_getMNIpos(patientID{id});
        for ii = 1: ns
            [label{ii},code{ii},roi_id(ii)]=fp_get_mni_anatomy(mni_pos(ii,:));
        end 
        u_roi_id = sort(unique(roi_id));
        nroi = numel(u_roi_id);
        
        for iroi = 2:nroi
            for ifq = 1:nfreq                
                ccoh = coh(:,roi_id == u_roi_id(iroi),roi_id == u_roi_id(iroi),ifq);
                rcoh(iit,:,iroi,iroi,ifq) = abs(imag(svd(ccoh))); % nit x npcs x nroi x nroi x nfq
                clear ccoh
                
            end
        end
        
    end
    
    clear rcoh_max
    rcoh_max = max(rcoh,2); %nroi x nroi x nfq
    COH = COH + rcoh_max;
    clear rcoh
end

%%%%here roi neighbouring 
freq_conn = zeros(nfreq,nfreq);
for ifreq = 1:nfreq-1
    freq_conn(ifreq,ifreq+1)=1;
    freq_conn(ifreq+1,ifreq)=1;
end
freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
kron_conn = kron(conn_s,freq_conn_s); %say nvox*nfreq = nkron
kron_conn = kron(conn_s,kron_conn);

threshold = prctile(reshape(COH,1,[]),99,9);
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


%shuffled clusters

onoff = COH>threshold;
big_clusters = zeros(nit-1,nfreq,ns,ns);

%find the clusters
for iit = 1: nit-1

    clear onoff_temp u ind NB
    onoff_temp = squeeze(onoff(iit,:,:)); %nfreq x ns
    u = onoff_temp(:); %should be the same indexing like in kron_conn now; nkron x 1

    ind = find(u==1); %remember indeces of super-threshold coherences
    NB = kron_conn;
    NB(u==0,:)=[]; %pass the neighbourhood structure only for the super-threshold voxels
    NB(:,u==0)=[];

    %components assigns every voxel to a cluster, even if this means that every voxel is its own cluster
    clear ci x clu
    [ci, x] = components(NB); %x is the histogram of the clusters
    clu = zeros(size(kron_conn,1),1);%refill with sub-threshold voxels
    clu(ind)= ci; %nkron x 1
    clu = reshape(clu,[nfreq ns ns]); %nfreq x ns

    if numel(x)>0
        big_clu_id = find(x==max(x));
        big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
        big_clusters(iit,:,:,:) = (clu == big_clu_id);
    end

end

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
