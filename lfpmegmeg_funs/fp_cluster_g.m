function [p, true_clu] = fp_cluster_g(minnbchan,fband, abs_imag,DIROUT)
%group statistics, cluster analysis with findclusters fun

fp_addpath 

if nargin>3
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

if strcmp(fband,'theta')
    frq_band = [4 8];
elseif strcmp(fband,'alpha')
    frq_band = [7 13];
elseif strcmp(fband,'beta')
    frq_band = [13 30];
elseif strcmp(fband,'gamma_low')
    frq_band = [30 60];
elseif strcmp(fband,'gamma_high')
    frq_band = [60 90];
else 
    warning('Choosing beta frequency band!')
    frq_band = [13 30];
    fband = 'beta';
end
    
fs = 300;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90) = [];
frq_id = find(frqs> frq_band(1) & frqs< frq_band(2));

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

avg_coh = squeeze(median(sum(COH(:,:,:,frq_id,:),1),4));
avg_coh = squeeze(reshape(avg_coh,[size(COH,2)*size(COH,3), size(COH,5)])); %all it x nvox
threshold = prctile(reshape(avg_coh,1,[]),99.9);
onoff = avg_coh>threshold;


%true cluster
clear clu total x big_clu_id
[clu, total] = findcluster(squeeze(onoff(1,:))',...
    match_conn, match_conn, minnbchan);
true_clu = clu';
true_total = total;
true_avg_coh = squeeze(avg_coh(1,:))';

onoff(1,:) =[]; %remove true coherence dimension
nit = size(onoff,1);


%shuffled clusters
big_clusters = zeros(nit,ns);
avg_coh = avg_coh(end-nit+1:end,:); %select shuffled clusters only

%find the clusters
for iit = 1: nit
    
    clear clu total x big_clu_id
    [clu, total] = findcluster(squeeze(onoff(iit,:))',...
        match_conn,  minnbchan);
    if total>0
        clear x
        x = hist(clu(:),0:total);
        big_clu_id = find(x(2:end)==max(x(2:end)));
        big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
        big_clusters(iit,:) = (clu == big_clu_id)';
    end
    
end

%compare not only cluster size but also magnitude of coherence within
%the relevant cluster
clear a
a = zeros(size(avg_coh)); %only shuffled clusters
a(big_clusters==1) = avg_coh(big_clusters==1);
shufCoh = squeeze(sum(a,2)); %cat across chunks

if true_total>0 %when at least one true cluster exists
    for iclus = 1:true_total
        clear trueCoh temp
        trueCoh = sum(true_avg_coh(true_clu==iclus)); %scalar
        p(iclus) = sum(shufCoh>trueCoh)/numel(shufCoh);
        if p(iclus) > 0.01
            break
        end
    end
    
elseif sum(shufCoh)== 0  %when no cluster was found it any iteration
    p= nan;
    
else %when only in shuffled conditions clusters were found
    clear trueCoh
    trueCoh = 0;
    p = sum(shufCoh>trueCoh)/numel(shufCoh);
end


outname = sprintf('%sp_cluster_g_%s_%s',DIROUT,fband,abs_imag);
save(outname,'p','true_clu','-v7.3')