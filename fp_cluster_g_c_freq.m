function [p, true_clu] = fp_cluster_g_c_freq(abs_imag, DIROUT)
%Group statistics, finds clusters with components fun.
%Clustering across space and frequencies.

fp_addpath

if nargin>2
    if ~exist(DIROUT); mkdir(DIROUT); end
else
    warning('Results will not be saved')
end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
nsubs = numel(patientID);

if isempty(abs_imag)
    abs_imag = 'abs';
end

[commonvox_pos, voxID] = fp_find_commonvox;
ns = size(commonvox_pos,1);

nchunk = 50;

for id = 1:nsubs
    
    %get neighbouring nodes and node positions
    clear conn mni_pos match_conn conn_s kron_conn sym_pos flip_id match_pos noEq
    mni_pos = fp_getMNIpos(patientID{id});
    conn = fp_find_neighbours(patientID{id}); %think about how to do that
    match_conn = conn(voxID{id},voxID{id});
    ns = size(match_conn,1);
    
    %get flip id and symmetric head
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    match_pos = sym_pos(voxID{id},:);
    [~,flip_id] = fp_flip_vol(match_pos);
    
    for ichunk = 1:nchunk
        %load coherences
        clear coh flip_coh abs_coh avg_coh
        load(sprintf('Coherences_Patient%s_chunk%d.mat',patientID{id},ichunk));
        
        [nit,nfreq,~,~] = size(coh);
        freq_conn = fp_get_freq_conn(nfreq);
        conn_s = sparse(match_conn); %nvox x nvox
        freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
        kron_conn = kron(conn_s,freq_conn_s); %say nvox*nfreq = nkron
        
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
        
        %median across lfp channels (already flipped) and across frequencies
        COH(id,ichunk,:,:,:) = squeeze(median(abs_coh,4));
    end
end

avg_coh = squeeze(sum(COH,1));
s_a_coh = reshape(avg_coh,[nchunk*nit,nfreq,ns]);
s_a_coh(1,:,:) = [];
threshold = prctile(reshape(s_a_coh,1,[]),99.9);
clear s_a_coh

%cat the chunks
avg_coh = squeeze(reshape(avg_coh,[size(COH,2)*size(COH,3),size(COH,4), size(COH,5)]));
onoff = avg_coh>threshold;


%true cluster
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

%save true cluster for later
clear true_clu true_total true_sizes
true_total = numel(x);
true_clu = fp_order_clusters(clu,true_total);
true_avg_coh = squeeze(avg_coh(1,:,:))';

onoff(1,:,:) =[]; %remove true coherence dimension
nit = size(onoff,1);


%shuffled clusters
shuf_clusters = zeros(nit,nfreq,ns);
avg_coh = avg_coh(end-nit+1:end,:,:);%select shuffled clusters only

%find the clusters
for iit = 1: nit
    
    clear onoff_temp u ind A
    onoff_temp = squeeze(onoff(iit,:,:)); %nfreq x ns
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
    total = numel(x);

    shuf_clusters(iit,:,:) = fp_order_clusters(clu,total);   
end


%compare not only cluster size but also magnitude of coherence within
%the relevant cluster

if true_total>0 %when at least one true cluster exists
    
   for iclus = 1:true_total
        
        clear a shufCoh
        a = zeros(size(avg_coh));
        a(shuf_clusters==iclus) = avg_coh(shuf_clusters==iclus);
        shufCoh = squeeze(sum(sum(a,2),3)); %cat across chunks
        
        clear trueCoh temp
        trueCoh = sum(sum(true_avg_coh(true_clu==iclus))); %scalar
        p(iclus) = sum(shufCoh>trueCoh)/numel(shufCoh);
    end
    
elseif sum(shuf_clusters)== 0  %when no cluster was found it any iteration
    p= nan;
    
else %when only in shuffled conditions clusters were found
    p = 1;
end


outname = sprintf('%sp_cluster_g_c_freq_%s',DIROUT,abs_imag);
save(outname,'p','threshold','true_clu','-v7.3')



