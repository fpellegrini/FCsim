function [p, true_clu] = fp_cluster_g_freq(minnbchan,abs_imag,DIROUT)
%group statistics, finds clusters across space and frequencies with the
%findclusters fun

fp_addpath

if nargin>2
    if ~exist(DIROUT); mkdir(DIROUT); end
else
    warning('Results will not be saved')
end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

if isempty(minnbchan)
    minnbchan = 2;
end
if isempty(abs_imag)
    abs_imag = 'abs';
end

[commonvox_pos, voxID] = fp_find_commonvox;

nchunk = 50;

for id = 1:numel(patientID)
    
    %get neighbouring nodes and node positions
    clear conn mni_pos match_conn sym_pos flip_id match_pos noEq
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
        nfreq = size(coh,2);
        
        %flip coherence
        coh(:,:,noEq,:) = [];
        match_coh = coh(:,:,voxID{id},:);
        flip_coh = match_coh;
        flip_coh(:,:,:,4:6) = match_coh(:,:,flip_id,4:6);
        
        %absolute value
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
threshold = prctile(reshape(avg_coh,1,[]),99.99);

%cat the chunks
avg_coh = squeeze(reshape(avg_coh,[size(COH,2)*size(COH,3),size(COH,4), size(COH,5)]));
onoff = avg_coh>threshold;

%true cluster
clear clu total x big_clu_id
[clu, true_total] = findcluster(squeeze(onoff(1,:,:))',...
    match_conn, match_conn, minnbchan);
true_clu = fp_order_clusters(clu',true_total);
true_avg_coh = squeeze(avg_coh(1,:,:))';

onoff(1,:,:) =[]; %remove true coherence dimension
nit = size(onoff,1);

%shuffled clusters
shuf_clusters = zeros(nit,ns,nfreq);
avg_coh = avg_coh(end-nit+1:end,:,:);%select shuffled clusters only

%find the clusters
for iit = 1: nit    
    clear clu total x big_clu_id
    [clu, total] = findcluster(squeeze(onoff(iit,:,:))',...
        match_conn,  minnbchan);
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
        if p(iclus) > 0.01
            break
        end
    end
    
elseif sum(shuf_clusters)== 0  %when no cluster was found it any iteration
    p= nan;
    
else %when only in shuffled conditions clusters were found
    p = 1;
end


outname = sprintf('%sp_cluster_g_freq_%s',DIROUT,abs_imag);
save(outname,'p','true_clu','-v7.3')