function threshold = fp_get_thresholds_group(fband, abs_imag, DIROUT)
%get threshold (from all chunks) for all subjects 

if nargin>2
    if ~exist(DIROUT); mkdir(DIROUT); end
end

if isempty(abs_imag)
    abs_imag = 'abs';
end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

if strcmp(fband,'theta')
    frq_band = [4 8];
elseif strcmp(fband,'alpha')
    frq_band = [7 13];
elseif strcmp(fband,'beta')
    frq_band = [13 30];
elseif strcmp(fband,'gamma_low')
    frq_band = [30 46];
elseif strcmp(fband,'gamma_high')
    frq_band = [60 90];
else 
    warning('Choosing beta frequency band!')
    frq_band = [13 30];
end

nchunk = 10;
fs = 300;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90) = [];
frq_id = find(frqs> frq_band(1) & frqs< frq_band(2));

coh_vals = [];

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
        avg_coh = squeeze(median(median(abs_coh(:,frq_id,:,:),4),2));
        coh_vals = cat(2,coh_vals, reshape(avg_coh,1,[]));
    end
end

clear threshold 
threshold = prctile(coh_vals,99);

outname = sprintf('%sthreshold_group_%s_%s',DIROUT,fband, abs_imag);
save(outname,'threshold','-v7.3')
    
    
    