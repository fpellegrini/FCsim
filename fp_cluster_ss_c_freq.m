function [p, TRUE_CLU] = fp_cluster_ss_c_freq(patientNumber,abs_imag,DIROUT)
%singlesub components with clustering across space and frequencies 

fp_addpath

if nargin>2
    if ~exist(DIROUT); mkdir(DIROUT); end
else
    warning('Results will not be saved')
end

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end
if isempty(abs_imag)
    abs_imag = 'abs';
end

nchunk = 10;

for id = 1:numel(patientID)
    
    %voxel con
    clear conn mni_pos noEq threshold
    conn = fp_find_neighbours(patientID{id});
    conn_s = sparse(conn); %nvox x nvox
    mni_pos = fp_getMNIpos(patientID{id});
    
    %get flip id and symmetric head
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    [~,flip_id] = fp_flip_vol(sym_pos);
    
    load(sprintf('%sthreshold_Patient%s_allfreq_%s.mat',DIROUT,patientID{id}, abs_imag));
    shufCoh = [];
    
    for ichunk = 1:nchunk
        
     %load coherences
        clear coh flip_coh abs_coh avg_coh
        load(sprintf('Coherences_Patient%s_chunk%d.mat',patientID{id},ichunk));
        
        coh(:,:,noEq,:) = [];
    [nit, nfreq, ns, nlfp] = size(coh);
        flip_coh = coh;
    flip_coh(:,:,:,4:6) = coh(:,:,flip_id,4:6);
       
    %freq conn and kron conn
    freq_conn = zeros(nfreq,nfreq);
    for ifreq = 1:nfreq-1
        freq_conn(ifreq,ifreq+1)=1;
        freq_conn(ifreq+1,ifreq)=1;        
    end   
    freq_conn_s = sparse(freq_conn);%nfreq x nfreq
    kron_conn = kron(conn_s,freq_conn_s); %say nvox*nfreq = nkron
    
    if strcmp(abs_imag,'abs')
        abs_coh= abs(flip_coh);
    elseif strcmp(abs_imag,'imag')
        abs_coh = abs(imag(flip_coh));
    else
        error('Method unknown!')
    end
    
    %median across lfp channels (already flipped) and across frequencies
    avg_coh = squeeze(median(abs_coh,4));
    
    clear onoff
    onoff = avg_coh>threshold;
    
    %true cluster
    if ichunk ==1 
        clear onoff_temp u ind A 
        onoff_temp = squeeze(onoff(1,:,:)); %nfreq x ns
        u = onoff_temp(:); %should be the same indexing like in kron_conn now; nkron x 1
        
        ind = find(u==1); %remember indeces of super-threshold coherences
        A = kron_conn;
        A(u==0,:)=[]; %pass the neighbourhood structure only for the super-threshold voxels
        A(:,u==0)=[];
        
        %components assigns every voxel to a cluster, even if this means that every voxel is its own cluster
        clear ci x clu
        [ci, x] = components(A); %x is the histogram of the clusters 
        clu = zeros(size(kron_conn,1),1);%refill with sub-threshold voxels
        clu(ind)= ci; %nkron x 1 
        clu = reshape(clu,[nfreq ns]); %nfreq x ns
         
        %save true cluster for later
            clear true_clu true_total true_sizes
            true_clu = clu;
            true_total = numel(x);
            true_sizes = x;
            true_avg_coh = squeeze(avg_coh(1,:,:))';
            
            onoff(1,:,:) =[]; %remove true coherence dimension
            nit = nit-1;      
    end
    
    %shuffled clusters 
     clear big_clusters
        big_clusters = zeros(nit,nfreq,ns);
        avg_coh = avg_coh(end-nit+1:end,:,:); %select shuffled clusters only
    
    for iit = 1: nit
        
        clear onoff_temp u ind A 
        onoff_temp = squeeze(onoff(iit,:,:)); %nfreq x ns
        u = onoff_temp(:); %should be the same indexing like in kron_conn now; nkron x 1
        
        ind = find(u==1); %remember indeces of super-threshold coherences
        A = kron_conn;
        A(u==0,:)=[]; %pass the neighbourhood structure only for the super-threshold voxels
        A(:,u==0)=[];
        
        %components assigns every voxel to a cluster, even if this means that every voxel is its own cluster
        clear ci x clu
        [ci, x] = components(A); %x is the histogram of the clusters 
        clu = zeros(size(kron_conn,1),1);%refill with sub-threshold voxels
        clu(ind)= ci; %nkron x 1 
        clu = reshape(clu,[nfreq ns]); %nfreq x ns      
            
      if numel(x)>0
            big_clu_id = find(x==max(x));
            big_clu_id=big_clu_id(1); %in case there are two clusters with the same size, take the first one
            big_clusters(iit,:,:) = (clu == big_clu_id);
      end        
    end
    
    %compare not only cluster size but also magnitude of coherence within
    %the relevant cluster 
    clear a
    a = zeros(size(avg_coh)); %only shuffled clusters
    a(big_clusters==1) = avg_coh(big_clusters==1);
        shufCoh = cat(1,shufCoh,squeeze(sum(sum(a,2),3))); %cat across chunks
    end
    
    
    if true_total>0 %when at least one true cluster exists  
        for iclus = 1:true_total
            clear trueCoh 
            trueCoh = sum(sum(true_avg_coh(true_clu==iclus)));
            p{id}(iclus) = sum(shufCoh>trueCoh)/numel(shufCoh);
        end
        TRUE_CLU{id} = true_clu;
        TRUE_sizes{id} = true_sizes;
        
    elseif sum(shufCoh)== 0  %when no cluster was found it any iteration
        p{id}= nan;
        
    else %when only in shuffled conditions clusters were found 
        clear trueCoh
        trueCoh = 0;
        p{id} = sum(shufCoh>trueCoh)/numel(shufCoh);
    end
    
end


if numel(patientID)==11 && exist('DIROUT')
    outname = sprintf('%sp_ss_c_allfreqs_%s',DIROUT,abs_imag);
    save(outname,'p','threshold','TRUE_CLU','TRUE_sizes','-v7.3')
end



