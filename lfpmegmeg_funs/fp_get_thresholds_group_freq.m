function threshold = fp_get_thresholds_group_freq(abs_imag, DIROUT)
%get threshold (from all chunks) for all subjects and all freqs

if nargin>2
    if ~exist(DIROUT); mkdir(DIROUT); end
end

if isempty(abs_imag)
    abs_imag = 'abs';
end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

nchunk = 10;

for id = 1:numel(patientID)  
    
    %get neighbouring nodes and node positions
    clear conn mni_pos noEq
    mni_pos = fp_getMNIpos(patientID{id});
    
    %get flip id and symmetric head
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    [~,flip_id] = fp_flip_vol(sym_pos); 
    
    %get flip id and symmetric head
    [~, noEq] = fp_symmetric_vol(mni_pos);
    
    for ichunk = 1:nchunk
        %load coherences
        clear coh flip_coh abs_coh avg_coh
        load(sprintf('Coherences_Patient%s_chunk%d.mat',patientID{id},ichunk));
        
        coh(:,:,noEq,:) = [];       
        flip_coh = coh;
        flip_coh(:,:,:,4:6) = coh(:,:,flip_id,4:6);
    
        if strcmp(abs_imag,'abs')
            abs_coh= abs(flip_coh);
        elseif strcmp(abs_imag,'imag')
            abs_coh = abs(imag(flip_coh));
        else
            error('Method unknown!')
        end
    
        %mean across lfp channels (already flipped) and across frequencies
        COH(id,ichunk,:,:,:) = squeeze(median(abs_coh,4));
        
    end
end

avg_coh = squeeze(sum(COH,1)); %sum over subjects
coh_vals = reshape(avg_coh,1,[]);

threshold = prctile(coh_vals,99);

outname = sprintf('%sthreshold_group_allfreq_%s',DIROUT, abs_imag);
save(outname,'threshold','-v7.3')
    
    
    