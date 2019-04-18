function p = fp_cluster_ss(patientNumber, fband, minnbchan,abs_imag)
%singlesubjects, finds clusters with the findclusters fun

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

fs = 300;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90) = [];
frq_id = find(frqs> frq_band(1) & frqs< frq_band(2));

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
    avg_coh = squeeze(median(median(abs_coh(:,frq_id,:,:),4),2));
    threshold(id) = prctile(reshape(avg_coh,1,[]),99);
    
    onoff = avg_coh>threshold(id);
    big_clusters = zeros(nit-1,ns);
    
    for iit = 1: nit
        
        clear clu total x big_clu_id
        [clu, total] = findcluster(squeeze(onoff(iit,:))',...
            conn, conn, minnbchan);
        
        if iit==1 %save true cluster for later
            true_clu = clu;
            true_total = total;
        elseif total>0
            clear x
            x = hist(clu,0:total);
            big_clu_id = find(x(2:end)==max(x(2:end)));
            big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
            big_clusters(iit-1,:) = clu == big_clu_id;
            
            %             if iit == 1
            %                 %plot biggest cluster
            %                 figure
            %                 c=sym_pos;
            %                 clear mask
            %     %             mask= onoff(iit,:);
            %                 mask = big_clusters(iit,:);
            %                 scatter3(c(:,1),c(:,2),c(:,3),5,[0.85 0.85 0.85])
            %                 hold on
            %                 scatter3(c(mask==1,1),c(mask==1,2),c(mask==1,3),20,[0.8 0.1,0.5],'filled')
            %                 colormap jet
            %
            %                 outname = sprintf('%s%s.png',DIRFIG,patientID{id});
            %                 print(outname,'-dpng');
            %                 close all
            %             end
            
            
        end
        
    end
    
    %compare not only cluster size but also magnitude of coherence within
    %the relevant cluster 
    b = avg_coh(2:end,:); %only shuffled clusters 
    a = zeros(size(b)); %only shuffled clusters
    a(big_clusters==1) = b(big_clusters==1);
    shufCoh = sum(a,2);
    
    
    if true_total>0 %when at least one true cluster exists  
        for iclus = 1:true_total
            clear trueCoh
            trueCoh = sum(avg_coh(1,true_clu==iclus));
            p{id}(iclus) = sum(shufCoh>trueCoh)/numel(shufCoh);
        end
        TRUE_CLU{id} = true_clu;
        
    elseif sum(shufCoh)== 0  %when no cluster was found it any iteration
        p{id}= nan;
        
    else %when only in shuffled conditions clusters were found 
        clear trueCoh
        trueCoh = 0;
        p{id} = sum(shufCoh>trueCoh)/numel(shufCoh);
    end
        
end

outname = sprintf('%sp_singlesub_%s',DIROUT,abs_imag);
save(outname,'p','threshold','TRUE_CLU','-v7.3')



