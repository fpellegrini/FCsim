function [p, true_clu] = fp_cluster_g_c_freq_test(abs_imag, DIROUT)
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
    abs_imag = 'imag';
end

[commonvox_pos, voxID] = fp_find_commonvox;
ns = size(commonvox_pos,1);
nchunk = 50;
%%
for id = 1:nsubs
    
    %get neighbouring nodes and node positions
    clear conn mni_pos match_conn conn_s kron_conn sym_pos flip_id match_pos noEq
    mni_pos = fp_getMNIpos(patientID{id});
    conn = fp_find_neighbours(patientID{id}); %think about how to do that
    match_conn = conn(voxID{id},voxID{id});
    
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

[nsub, ~,~,nfreq,ns] = size(COH);
nCOH = squeeze(reshape(COH,[nsub,size(COH,2)*size(COH,3),nfreq,ns]));
tCoh = squeeze(nCOH(:,1,:,:));
sCoh = nCOH(:,2:end,:,:);
nit = size(sCoh,2);

%% true 

clear o cCoh csCoh dbCoh ps hs testval

%debias
tic
o = squeeze(mean(sCoh,2)); 
dbCoh = tCoh-o; 

%sign-rank test
for ifreq = 1:nfreq
    for ivox = 1:ns
%         [hn(ivox,ifreq), pn(ivox,ifreq)] = kstest(dbCoh(:,ifreq,ivox)); %not one is n.d.
%         [ht(ivox,ifreq), pt(ivox,ifreq),stats] = ttest(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',0.001);
%         testval = stats(1);
        [ps(ivox,ifreq), hs(ivox,ifreq),stats] = signrank(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',0.001);
        testval(ivox,ifreq) = stats.signedrank;
    end
end
toc

onoff = hs;
clear u ind A

u = onoff(:); %should be the same indexing like in kron_conn now; nkron x 1
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

true_total = numel(x);
true_clu = fp_order_clusters(clu,true_total);
true_testval = testval;
true_p = ps;


%% shuffled

for iit = 1:nit 
    tic
    clear o cCoh csCoh dbCoh ps hs stats testval onoff 
    
    %select current shuffled coh
    cCoh = squeeze(sCoh(:,iit,:,:)); 
    if iit == 1
        csCoh = sCoh(:,2:end,:,:); 
    elseif iit == nit
        csCoh = sCoh(:,1:end-1,:,:);
    else 
        csCoh = sCoh(:,[1:iit-1 iit+1:end],:,:);
    end
    
    %debias
    o = squeeze(mean(csCoh,2)); 
    dbCoh = cCoh-o; 

    %sign-rank test
    for ifreq = 1:nfreq
        for ivox = 1:ns
    %         [hn(ivox,ifreq), pn(ivox,ifreq)] = kstest(dbCoh(:,ifreq,ivox)); %not one is n.d.
    %         [ht(ivox,ifreq), pt(ivox,ifreq),stats] = ttest(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',0.001);
    %         testval = stats(1);
            [ps(ivox,ifreq), hs(ivox,ifreq),stats] = signrank(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',0.001);
            testval(ivox,ifreq) = stats.signedrank;
        end
    end

    onoff = hs;
    clear u ind A
    u = onoff(:); %should be the same indexing like in kron_conn now; nkron x 1
    
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
    shuf_clu(iit,:,:) = fp_order_clusters(clu,total);
    shuf_total(iit) = total;
    shuf_testval(iit,:,:) = testval;
    shuf_p(iit,:,:) = ps;
    
    toc
end

%% 

if true_total>0 %when at least one true cluster exists    
    
    for iclus = 1:true_total
        
        clear a shuf_clu_val
        a = zeros(size(shuf_testval));
%         a(shuf_clu==iclus) = shuf_testval(shuf_clu==iclus); %compare first clu with first, second with second etc
        a(shuf_clu==1) = shuf_testval(shuf_clu==1); %take always the biggest cluster 
        shuf_clu_val = squeeze(sum(sum(a,2),3));
        
        clear trueCoh temp
        true_clu_val = sum(sum(true_testval(true_clu==iclus))); %scalar
        p(iclus) = sum(shuf_clu_val>true_clu_val)/numel(shuf_clu_val);
        if p(iclus) < 0.01
            break
        end
    end
    
elseif sum(shuf_clusters)== 0  %when no cluster was found it any iteration
    p= nan;
    
else %when only in shuffled conditions clusters were found
    p = 1;
end
%%
outname = sprintf('%sp_cluster_g_c_freq_%s_test',DIROUT,abs_imag);
save(outname,'p','true_clu','true_p','-v7.3')



