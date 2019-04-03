function p = fp_cluster_group(threshold_coh, minnbchan,abs_imag)

cd ~/Dropbox/MEG_Project/Data

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

if isempty(threshold_coh)
    threshold_coh = 0.085;
end
if isempty(minnbchan)
    minnbchan = 2;
end 
if isempty(abs_imag)
    abs_imag = 'abs';
end 

nsubs = numel(patientID);
nit = 51;
nfreq = 46;
ns = 1917; %this is the issue
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
    
    %get flip id and symmetric head
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    [~,flip_id] = fp_flip_vol(sym_pos);
    
    %flip coherence
    coh(:,:,noEq,:) = [];
    flip_coh = coh;
    flip_coh(:,:,:,4:6) = coh(:,:,flip_id,4:6);
    
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
onoff = avg_coh>threshold_coh;

big_clusters = zeros(nit,nfreq,ns);

for iit = 1: nit

    clear cluster total x big_clu_id
    [cluster, total] = findcluster(squeeze(onoff(iit,:,:))',conn,conn, minnbchan); %try components fun

    if total>0
        x = hist(cluster(:),0:total);
        big_clu_id = find(max(x));
        big_clusters(iit,:,:) = cluster == big_clu_id;
    end

end

a = zeros(size(avg_coh));
a(big_clusters==1)=avg_coh(big_clusters==1);
clustercoh = sum(sum(a,2),3);

shufCoh = clustercoh(2:end);
trueCoh = clustercoh(1);
p = sum(shufCoh>trueCoh)/numel(shufCoh);



