function [p, true_clu] = fp_cluster_g_c_freq(abs_imag)
%Group statistics, finds clusters with components fun.
%Clustering across space and frequencies.

DIROUT = '~/Dropbox/Data_MEG_Project/';
if ~exist(DIROUT); mkdir(DIROUT); end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

if isempty(abs_imag)
    abs_imag = 'abs';
end

[commonvox_pos, voxID] = fp_find_commonvox;

nsubs = numel(patientID);
nit = 51;
nfreq = 46;
ns = size(commonvox_pos,1);

COH = nan(nsubs,nit,nfreq,ns);

for id = 1:nsubs
    
    %load coherences
    clear coh
    load(sprintf('Coherences_Patient%s.mat',patientID{id}));
    
    %get neighbouring nodes and node positions
    clear conn mni_pos match_conn conn_s kron_conn sym_pos flip_id match_pos noEq
    conn = fp_find_neighbours(patientID{id});
    match_conn = conn(voxID{id},voxID{id});
    freq_conn = zeros(nfreq,nfreq);
    for ifreq = 1:nfreq-1
        freq_conn(ifreq,ifreq+1)=1;
        freq_conn(ifreq+1,ifreq)=1;
    end
    conn_s = sparse(match_conn); %nvox x nvox
    freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
    kron_conn = kron(conn_s,freq_conn_s); %say nvox*nfreq = nkron
    
    %get flip id and symmetric head
    mni_pos = fp_getMNIpos(patientID{id});
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    match_pos = sym_pos(voxID{id},:);
    [~,flip_id] = fp_flip_vol(match_pos);
    
    coh(:,:,noEq,:) = [];
    match_coh = coh(:,:,voxID{id},:);
    flip_coh = match_coh;
    flip_coh(:,:,:,4:6) = coh(:,:,flip_id,4:6);
    
    if strcmp(abs_imag,'abs')
        abs_coh= abs(flip_coh);
    elseif strcmp(abs_imag,'imag')
        abs_coh = abs(imag(flip_coh));
    else
        error('Method unknown!')
    end
    
    %mean across lfp channels
    COH(id,:,:,:) = squeeze(mean(abs_coh,4));
end

clear avg_coh threshold onoff big_clusters
avg_coh = squeeze(sum(COH,1)); %sum across subjects
threshold(id) = prctile(reshape(avg_coh,1,[]),99);
onoff = avg_coh>threshold(id);
big_clusters = zeros(nit-1,nfreq,ns);

for iit = 1: nit
    
    clear onoff_temp u ind A
    onoff_temp = squeeze(onoff(1,:,:)); %nfreq x ns
    u = onoff_temp(:); %should be the same indexing like in kron_conn now; nkron x 1
    
    ind = find(u==1); %remember indeces of super-threshold coherences
    A = kron_conn;
    A(u==0,:)=[]; %pass the neighbourhood structure only for the super-threshold voxels
    A(:,u==0)=[];
    
    %components assigns every voxel to a cluster, even if this means that every voxel is its own cluster
    clear ci x clu
    [ci, x] = components(A); %x is the histogram of the clusters
    clu = zeros(size(kron_conn,1),1);%refill with sub-threshold voxels
    clu(ind)= ci; %nkron x 1
    clu = reshape(clu,[nfreq ns]); %nfreq x ns
    
    if iit==1 %save true cluster for later
        clear true_clu true_total true_sizes
        true_clu = clu;
        true_total = numel(x);
        true_sizes = x;
        
    elseif numel(x)>0
        big_clu_id = find(x==max(x));
        big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
        big_clusters(iit-1,:,:) = (clu == big_clu_id);
    end
    
end

%compare not only cluster size but also magnitude of coherence within
%the relevant cluster
clear b a shufCoh
b = avg_coh(2:end,:,:); %only shuffled clusters
a = zeros(size(b)); %only shuffled clusters
a(big_clusters==1) = b(big_clusters==1);
shufCoh = sum(sum(a,2),3); %nit x 1


if true_total>0 %when at least one true cluster exists
    for iclus = 1:true_total
        clear trueCoh temp
        temp = squeeze(avg_coh(1,:,:)); %select only true coherence
        trueCoh = sum(sum(temp(true_clu==iclus))); %scalar
        p(iclus) = sum(shufCoh>trueCoh)/numel(shufCoh);
    end
    
elseif sum(shufCoh)== 0  %when no cluster was found it any iteration
    p= nan;
    
else %when only in shuffled conditions clusters were found
    clear trueCoh
    trueCoh = 0;
    p = sum(shufCoh>trueCoh)/numel(shufCoh);
end


outname = sprintf('%sp_group_c_allfreqs_%s',DIROUT,abs_imag);
save(outname,'p','threshold','true_clu','true_sizes','-v7.3')



