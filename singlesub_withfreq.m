function p = singlesub_withfreq(patientNumber, minnbchan,abs_imag)

cd ~/Dropbox/Data_MEG_Project/

DIROUT = '~/Dropbox/Data_MEG_Project/';
if ~exist(DIROUT); mkdir(DIROUT); end

DIRFIG = sprintf('~/Dropbox/Data_MEG_Project/figures/cluster_singlesub/%s/',abs_imag);
if ~exist(DIRFIG); mkdir(DIRFIG); end

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

if isempty(minnbchan)
    minnbchan = 2;
end
if isempty(abs_imag)
    abs_imag = 'abs';
end

fs = 300;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90) = [];

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
    
    if strcmp(abs_imag,'abs')
        abs_coh= abs(flip_coh);
    elseif strcmp(abs_imag,'imag')
        abs_coh = abs(imag(flip_coh));
    else
        error('Method unknown!')
    end
    
    %mean across lfp channels (already flipped) and across frequencies
    avg_coh = squeeze(median(abs_coh(:,:,:,:),4));
    avg_coh = permute(avg_coh,[3 2 1]);
%     for i = 1:46
%         avg_coh(:,i,:) = zscore(avg_coh(:,i,:));
%     end
    threshold(id) = prctile(reshape(avg_coh,1,[]),99);
    
    onoff = avg_coh>threshold(id);
    big_clusters = zeros(nit,ns,nfreq);
    
    for iit = 1: nit
        
        clear cluster total x big_clu_id
        [clu, total] = findcluster(squeeze(onoff(:,:,iit)),...
            conn, conn, minnbchan);
        
        if total>0
            clear x
            x = hist(clu(:),0:total);
            big_clu_id = find(x(2:end)==max(x(2:end)));
            big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
            big_clusters(iit,:,:) = clu == big_clu_id;
            
         
            
        end
        
    end
    
    a = zeros(size(avg_coh));
    a(big_clusters==1)=avg_coh(big_clusters==1);
    clustercoh = sum(a,2);
    
    if sum(clustercoh)> 0
        shufCoh = clustercoh(2:end);
        trueCoh = clustercoh(1);
        p(id) = sum(shufCoh>trueCoh)/numel(shufCoh);
    else %when no cluster was found in any iteration 
        p(id)= nan;
    end
    
end





