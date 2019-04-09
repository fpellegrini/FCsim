function p = fp_cluster_group(minnbchan,abs_imag)

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
fs = 300;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90) = [];

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

avg_coh = squeeze(median(sum(COH(:,:,13:30,:),1),3));
threshold = prctile(reshape(avg_coh,1,[]),99);
onoff = avg_coh>threshold;

% big_clusters = zeros(nit,nfreq,ns);
big_clusters = zeros(nit,ns);

for iit = 1: nit

    clear cluster total x big_clu_id
    [clu, total] = findcluster(squeeze(onoff(iit,:))',match_conn,match_conn, minnbchan);

    if total>0
        x = hist(clu(:),0:total);
        big_clu_id = find(x(2:end)==max(x(2:end)));
        big_clusters(iit,:,:) = clu == big_clu_id;
        
%         figure
%         c=commonvox_pos;
%         clear mask
%         mask = big_clusters(iit,:);
%         scatter3(c(:,1),c(:,2),c(:,3),5,[0.85 0.85 0.85])
%         hold on
%         scatter3(c(mask==1,1),c(mask==1,2),c(mask==1,3),20,[0.8 0.1,0.5],'filled')
%         colormap jet
    end

end

a = zeros(size(avg_coh));
a(big_clusters==1)=avg_coh(big_clusters==1);
clustercoh = sum(sum(a,2),3);

shufCoh = clustercoh(2:end);
trueCoh = clustercoh(1);
p = sum(shufCoh>trueCoh)/numel(shufCoh);

outname = sprintf('%sp_group_%s',DIROUT,abs_imag);
save(outname,'p','-v7.3')