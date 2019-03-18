function fp_cluster_coh_beta(patientNumber)

cd ~/Dropbox/MEG_Project/Data

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

minnbchan = 1;

fs = D.fsample;
fres = fs;
frqs = sfreqs(fres, fs);
frq_inds = find(frqs > 13 & frqs < 30);

for id = 1:numel(patientID)
    
    load(sprintf('Coherences_Patient%s.mat',patientID{id}));
    
    coh = mean(coh(:,frq_inds,:,:),2);
    [~, nfreq, ns, nlfp] = size(coh);

    for ifreq = 1:nfreq
        for is = 1:ns
            for ilfp = 1:nlfp
                shufCoh = squeeze(abs(coh(2:end,ifreq,is,ilfp)));
                trueCoh = squeeze(abs(coh(1,ifreq,is,ilfp)));
                p(ifreq,is,ilfp) = sum(shufCoh>trueCoh)/numel(shufCoh);
                clear shufCoh trueCoh
            end 
        end
    end

    onoff = p < 0.95;
    
    onoff1 = squeeze(onoff(:,:,1))';
    conn = fp_find_neighbours(patientID{id});
    
    [cluster, total] = findcluster(onoff1, conn,conn, minnbchan);
    
    a = find(cluster==1);
    scatter3(c(:,1),c(:,2),c(:,3))
    hold on
    scatter3(c(a,1),c(a,2),c(a,3),'r+')
    
    
    
end

