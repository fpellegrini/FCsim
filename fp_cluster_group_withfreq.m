function [p, true_clu] = fp_cluster_group_withfreq(minnbchan,abs_imag)

cd ~/Dropbox/Data_MEG_Project/

DIROUT = '~/Dropbox/Data_MEG_Project/';
if ~exist(DIROUT); mkdir(DIROUT); end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

if isempty(minnbchan)
    minnbchan = 2;
end
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
    load(sprintf('Coherences_Patient%s.mat',patientID{id}));
    mni_pos = fp_getMNIpos(patientID{id});
    conn = fp_find_neighbours(patientID{id}); %think about how to do that
    match_conn = conn(voxID{id},voxID{id});
    
    %get flip id and symmetric head
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    match_pos = sym_pos(voxID{id},:);
    [~,flip_id] = fp_flip_vol(match_pos);
    
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
    
    %sum across lfp channels
    COH(id,:,:,:) = squeeze(mean(abs_coh,4));
    
end

avg_coh = squeeze(sum(COH,1));
threshold = prctile(reshape(avg_coh,1,[]),99);
onoff = avg_coh>threshold;

big_clusters = zeros(nit-1,nfreq,ns);
% big_clusters = zeros(nit,ns);

%find the clusters 
for iit = 1: nit
    
    clear clu total x big_clu_id
    [clu, total] = findcluster(squeeze(onoff(iit,:,:))',...
        match_conn,  minnbchan);
    if iit==1 %save true cluster for later
        true_clu = clu';
        true_total = total;
    elseif total>0
        clear x
        x = hist(clu(:),0:total);
        big_clu_id = find(x(2:end)==max(x(2:end)));
        big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
        big_clusters(iit-1,:,:) = (clu == big_clu_id)';
    end
    
end

%compare not only cluster size but also magnitude of coherence within
%the relevant cluster
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


outname = sprintf('%sp_group_withfreq_%s',DIROUT,abs_imag);
save(outname,'p','true_clu','-v7.3')