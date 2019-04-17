function p = singlesub_components(patientNumber, frq_band, minnbchan,abs_imag)
%doesnt work yet 
cd ~/Dropbox/Data_MEG_Project/

DIROUT = '~/Dropbox/Data_MEG_Project/';
if ~exist(DIROUT); mkdir(DIROUT); end

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

if isempty(frq_band)
    frq_band = [13 30];
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
frq_id = find(frqs> frq_band(1) & frqs< frq_band(2));

for id = 5 %1:numel(patientID)
    
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
    avg_coh = squeeze(median(median(abs_coh(:,frq_id,:,:),4),2));
    threshold(id) = prctile(reshape(avg_coh,1,[]),99);
    
    clear onoff
    onoff = avg_coh>threshold(id);
    big_clusters = zeros(nit,ns);
    
    for iit = 1: nit
        
            
        ind = find(onoff==1);
        A = conn;
        A(onoff(iit,:)'==0,:)=[];
        A(:,onoff(iit,:)'==0)=[];
        P = sparse(A);
        [ci sizes] = components(P);
        
%         clear cluster total x big_clu_id
%         [clu, total] = findcluster(squeeze(onoff(iit,:))',...
%             conn, conn, minnbchan);
%         
        if total>0
            clear x
            x = hist(clu,0:total);
            big_clu_id = find(max(x(2:end)));
            big_clusters(iit,:) = clu == big_clu_id;
            
            keyboard
            figure
            load('BF_Patient10.mat')
            clear mask
            mask = onoff(iit,:);
            c = sources.grid.pos;
            scatter3(c(:,1),c(:,2),c(:,3),5,[0.85 0.85 0.85])
            hold on
            scatter3(c(mask==1,1),c(mask==1,2),c(mask==1,3),20,[0.8 0.1,0.5],'filled')
            colormap jet
            
        end
        
    end
    
    a = zeros(size(avg_coh));
    a(big_clusters==1)=avg_coh(big_clusters==1);
    clustercoh = sum(a,2);
    
    shufCoh = clustercoh(2:end);
    trueCoh = clustercoh(1);
    p(id) = sum(shufCoh>trueCoh)/numel(shufCoh);
    
end

outname = sprintf('%sp_singlesub_%s1',DIROUT,abs_imag);
save(outname,'p','threshold','-v7.3')



