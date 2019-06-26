function fp_megmeg_pipeline(patientNumber)

% fp_addpath

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
    
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    X = D(:,:,:);
    X = X./10^(log10(range(X(:)))-2);
    
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
    
    CS = CS(1:(end-nlfp),1:(end-nlfp),:); %throw away lfp channels
    
    %construct filters
    
    %leadfield
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    L1 = inverse.MEG.L;
    ns = numel(L1);
    for is=1:ns
        L(:,is,:)= L1{is};
    end
    L = L.* (10^(-log10(range(L(:)))));
    
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
    
    %calculate coherence for permutations
    
    for iit = 1:nit
        
        %cross spectrum
        clear CS
        id_trials_1 = 1:n_trials;
        rng('shuffle')
        id_trials_2 = randperm(n_trials);
        CS = fp_tsdata_to_cpsd(X,fres,'MT',id_meg_chan, id_meg_chan, id_trials_1, id_trials_2);
        
        
        clear mni_pos label code roi_id u_roi_id
        mni_pos = fp_getMNIpos(patientID{id});
        for ii = 1: ns
            [label{ii},code{ii},roi_id(ii)]=fp_get_mni_anatomy(mni_pos(ii,:));
        end
        u_roi_id = sort(unique(roi_id));
        nroi = numel(u_roi_id);
        
        csroi = nan(nroi,nfreq,5,5);
        %project cross spectrum to voxel space and get power and coherence
        
                
        clear v5
        for ifq = 1:nfreq
            
            clear Aroi A_ CSv pv CSn
            Aroi = squeeze(A(:,:,roi_id~=0,ifq));
            A_ = reshape(Aroi, [nmeg, 3*size(Aroi,3)]);
            CSv = A_' * CS(:,:,ifq) * A_;
            pv = real(diag(CSv)); %3*nvox x 1
            CSn = CSv ./ sqrt(pv * pv');
            
            %region pca 
            is = 1;
            for iroi = 2:nroi
                clear v cCS cns
                
                cns = sum(roi_id == u_roi_id(iroi))*3;
                cCS = CSn(is:is+cns-1,is:is+cns-1);
                [v, ~, ~] = eig(real(cCS));
                
                if size(v,1)>5
                    v5(is:is+size(v,1)-1,:,ifq) = v(:,1:5); %5 * nregionvoxels
                else
                    v5(is:is+size(v,1)-1,1:size(v,2),ifq) = v;
                end
                
                is = is+cns;
            end
            
            %apply the filters to the cs 
            kr=1;
            for kroi = 2:nroi
                
                vk = v5(kr: kr+ sum(roi_id == u_roi_id(kroi))-1,:,ifq);
                jr=1;
                for jroi = 2: nroi
                    clear cCS
                    cCS = CSn(roi_id == u_roi_id(kroi),roi_id == u_roi_id(jroi));
                    vj = v5(jr: jr+ sum(roi_id == u_roi_id(jroi))-1,:,ifq);
                    cseig(kroi,jroi,:,:) = vk' * cCS * vj;
                    
                    jr=jr+sum(roi_id == u_roi_id(jroi));
                end
                
                kr=kr+sum(roi_id == u_roi_id(kroi));
            end
            
            %divide by power 
            for ifc=1:5
                for jfc =1:5
                    clear proi 
                    proi = squeeze(real(diag(cseig(:,:,ifc,jfc))));
                    csroi(:,:,ifc,jfc,ifq)= cseig(:,:,ifc,jfc)./sqrt(proi' * proi);
                end
            end
            
            
        end
        
        
        
    
        
        %%%% + cca , also dann coherence nit x nroi x nroi x nfreq (jeweils
        %%%%schon ueber ids aufsummieren!)
        
        
        COH(iit,:,:,:) = COH(iit,:,:,:) + coh;
        
    end
end


%neighbourhood
load('roi_conn.mat')
roiconn_s = sparse(roi_conn);
freq_conn = fp_get_freq_conn(nfreq);
freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
kron_conn = kron(roiconn_s,freq_conn_s); %nroi*nfreq = nkron
kron_conn = kron(roiconn_s,kron_conn);

threshold = prctile(reshape(COH,1,[]),99,9);


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
