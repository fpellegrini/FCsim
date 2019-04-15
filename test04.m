patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

ALL = [];

for id = 1:numel(patientID)
    
    %load coherences
    load(sprintf('Coherences_Patient%s.mat',patientID{id}));
    
    %get neighbouring nodes and node positions
    clear conn mni_pos noEq
    conn = fp_find_neighbours(patientID{id});
    mni_pos = fp_getMNIpos(patientID{id});
    
    %get flip id and symmetric head
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    [~,flip_id] = fp_flip_vol(sym_pos);
    
    coh(:,:,noEq,:) = [];
    [nit, nfreq, ns, nlfp] = size(coh);
    flip_coh = coh;
    flip_coh(:,:,:,4:6) = coh(:,:,flip_id,4:6);
    
    
    abs_coh= abs(flip_coh);
    
    
    %mean across lfp channels (already flipped) and across frequencies
    avg_coh = squeeze(median(abs_coh(:,:,:,:),4));
    avg_coh = permute(avg_coh,[3 2 1]);
    
%     a = avg_coh(:,:,1);
%     figure
%     imagesc(a)
%     caxis([0 0.3])
%     colorbar
    
    threshold(id) = prctile(reshape(avg_coh,1,[]),95);
    
    ALL = cat(2, ALL,reshape(avg_coh,1,[]));
    clearvars -except patientID id threshold ALL
end 

thresh_all = prctile(ALL,95);