function fp_cluster_coh(patientNumber)

cd ~/Dropbox/MEG_Project/Data

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

for id = 1:numel(patientID)
    
    load(sprintf('Coherences_Patient%s.mat',patientID{id}));
    
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

    clust = p < 0.95;
    
    
    
end

