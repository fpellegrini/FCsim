function p = singlesub_components_withfreq(patientNumber,abs_imag)
%doesnt work yet 
cd ~/Dropbox/Data_MEG_Project/

DIROUT = '~/Dropbox/Data_MEG_Project/';
if ~exist(DIROUT); mkdir(DIROUT); end

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end
if isempty(abs_imag)
    abs_imag = 'abs';
end


for id = 1 %1:numel(patientID)
    
    %load coherences
    clear coh
    load(sprintf('Coherences_Patient%s_chunk1.mat',patientID{id}));
    nfreq = size(coh,2);
    
    %get neighbouring nodes and node positions
    clear conn mni_pos noEq
    conn = fp_find_neighbours(patientID{id});  
    freq_conn = zeros(nfreq,nfreq);
    for ifreq = 1:nfreq-1
        freq_conn(ifreq,ifreq+1)=1;
        freq_conn(ifreq+1,ifreq)=1;        
    end
    conn_s = sparse(conn); %nvox x nvox
    freq_conn_s = sparse(freq_conn);%nfrerq x nfreq
    kron_conn = kron(conn_s,freq_conn_s); %say nvox*nfreq = nkron 
    
    %get flip id and symmetric head
    mni_pos = fp_getMNIpos(patientID{id});
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
    avg_coh = squeeze(median(abs_coh,4));
    threshold(id) = prctile(reshape(avg_coh,1,[]),99);
    
    clear onoff
    onoff = avg_coh>threshold(id);
    big_clusters = zeros(nit-1, nfreq,ns);
    
    for iit = 1: nit
        
        onoff_temp = squeeze(onoff(1,:,:));
        u = onoff_temp(:); %should be the same indexing like in kron_conn now; nkron x 1
        
        ind = find(u==1); %remember indeces of super-threshold coherences
        A = kron_conn;
        A(u==0,:)=[]; %pass the neighbourhood structure only for the super-threshold voxels
        A(:,u==0)=[];
        
        %components assigns every voxel to a cluster, even if this means that every voxel is its own cluster
        [ci, x] = components(A); %x is the histogram of the clusters 
        clu = zeros(size(kron_conn,1),1);%refill with sub-threshold voxels
        clu(ind)= ci; %nkron x 1 
         
        if iit==1 %save true cluster for later
            true_clu = clu;
            true_total = numel(x);
            true_siz = x;
            
        elseif numel(x)>0
            big_clu_id = find(x==max(x));
            big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
            big_clusters(iit-1,:,:) = (clu == big_clu_id)';
        end
        
    end
    
    %%modify from here 
    
    %compare not only cluster size but also magnitude of coherence within
    %the relevant cluster 
    b = avg_coh(2:end,:,:); %only shuffled clusters 
    a = zeros(size(b)); %only shuffled clusters
    a(big_clusters==1) = b(big_clusters==1);
    shufCoh = sum(sum(a,2),3); %nit x 1
    
    
    if true_total>0 %when at least one true cluster exists  
        for iclus = 1:true_total
            clear trueCoh temp
            temp = squeeze(avg_coh(1,:,:)); %select only true coherence
            trueCoh = sum(sum(temp(true_clu==iclus))); %scalar
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

outname = sprintf('%sp_singlesub_c_allfreqs_%s',DIROUT,abs_imag);
save(outname,'p','threshold','TRUE_CLU','-v7.3')



