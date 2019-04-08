function fp_cluster_coh_beta(patientNumber)

cd ~/Dropbox/Data_MEG_Project/

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

minnbchan = 1;
lfpchan = 2;

fs = D.fsample;
fres = 75;
frqs = sfreqs(fres, fs);
frq_inds = find(frqs > 12 & frqs < 31);

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

    onoff = p < 0.01;
    
    for lfpchan = 1:6
        figure
        onoff1 = squeeze(onoff(:,:,lfpchan))';
        conn = fp_find_neighbours(patientID{id});

        [cluster, total] = findcluster(onoff1, conn,conn, minnbchan);

        if total>0
            x = hist(cluster,0:total);
            cluSize = max(x(2:end));
            mask = cluster==(find(x==max(x(2:end)))-1);
        end

        c = sources.grid.pos;
        scatter3(c(:,1),c(:,2),c(:,3),5,[0.85 0.85 0.85])
        hold on
        col = -log(squeeze(p(:,:,lfpchan)));
        scatter3(c(mask,1),c(mask,2),c(mask,3),20,col(mask)','filled')
        colormap jet
    end 
    
    
    
end

