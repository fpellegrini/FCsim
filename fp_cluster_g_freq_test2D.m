function [p, true_clu] = fp_cluster_g_freq_test2D(minnbchan,abs_imag,DIROUT)
%group statistics, finds clusters across space and frequencies with the
%findclusters fun

fp_addpath

if nargin>2
    if ~exist(DIROUT); mkdir(DIROUT); end
else
    warning('Results will not be saved')
end

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
nsubs = numel(patientID);

if isempty(minnbchan)
    minnbchan = 2;
end
if isempty(abs_imag)
    abs_imag = 'imag';
end

[~, voxID] = fp_find_commonvox;
nchunk = 50;
alpha = 0.001;
%%
for id = 1:nsubs
    fprintf('Working on subject %s',patientID{id})
    
    %get neighbouring nodes and node positions
    clear mni_pos sym_pos flip_id match_pos noEq
    
    %get flip id and symmetric head
    mni_pos = fp_getMNIpos(patientID{id});
    [sym_pos, noEq] = fp_symmetric_vol(mni_pos);
    match_pos = sym_pos(voxID{id},:);
    [~,flip_id] = fp_flip_vol(match_pos);
    
    for ichunk = 1:nchunk
        %load coherences
        clear coh flip_coh abs_coh avg_coh
        load(sprintf('Coherences_e2D_Patient%s_chunk%d.mat',patientID{id},ichunk));
        
        %flip coherence
        coh(:,:,noEq,:) = [];
        match_coh = coh(:,:,voxID{id},:);
        flip_coh = match_coh;
        flip_coh(:,:,:,4:6) = match_coh(:,:,flip_id,4:6);
        
        %absolute value
        if strcmp(abs_imag,'abs')
            r_coh= log10(abs(flip_coh));
        elseif strcmp(abs_imag,'imag')
            r_coh = log10(abs(imag(flip_coh)));
            r_coh(:,1,:,:) = []; %delete inf at freq=1
        else
            error('Method unknown!')
        end
        
        %median across lfp channels (already flipped) and across frequencies
        COH(id,ichunk,:,:,:) = squeeze(median(r_coh,4));       
    end
end

[nsub, nchunk,niit,nfreq,ns] = size(COH);

conn = fp_find_neighbours(patientID{id});
match_conn = conn(voxID{id},voxID{id});

nCOH = reshape(COH,[nsub,nchunk*niit,nfreq,ns]);
tCoh = squeeze(nCOH(:,1,:,:));
sCoh = nCOH(:,2:end,:,:);
nit = size(sCoh,2);

%% true
 
clear o cCoh csCoh dbCoh ps hs testval

%debias
tic
o = squeeze(mean(sCoh,2)); 
dbCoh = tCoh-o; 

%sign-rank test
for ifreq = 1:nfreq
    for ivox = 1:ns
%         [hn(ivox,ifreq), pn(ivox,ifreq)] = kstest(dbCoh(:,ifreq,ivox)); %not one is n.d.
%         [ht(ivox,ifreq), pt(ivox,ifreq),stats] = ttest(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',0.001);
%         testval = stats(1);
        [ps(ivox,ifreq), hs(ivox,ifreq),stats] = signrank(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',alpha);
        testval(ivox,ifreq) = stats.signedrank;
    end
end
toc

onoff = hs;
clear clu total
[clu, total] = findcluster(onoff,...
    match_conn, match_conn, minnbchan);
true_clu = fp_order_clusters(clu',total);
true_total = total;
true_testval = testval;
true_p = ps;

%% shuffled


for iit = 1:nit 
    tic
    clear o cCoh csCoh dbCoh ps hs stats testval onoff 
    
    %select current shuffled coh
    cCoh = squeeze(sCoh(:,iit,:,:)); 
    if iit == 1
        csCoh = sCoh(:,2:end,:,:); 
    elseif iit == nit
        csCoh = sCoh(:,1:end-1,:,:);
    else 
        csCoh = sCoh(:,[1:iit-1 iit+1:end],:,:);
    end
    
    %debias
    o = squeeze(mean(csCoh,2)); 
    dbCoh = cCoh-o; 

    %sign-rank test
    for ifreq = 1:nfreq
        for ivox = 1:ns
    %         [hn(ivox,ifreq), pn(ivox,ifreq)] = kstest(dbCoh(:,ifreq,ivox)); %not one is n.d.
    %         [ht(ivox,ifreq), pt(ivox,ifreq),stats] = ttest(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',0.001);
    %         testval = stats(1);
            [ps(ivox,ifreq), hs(ivox,ifreq),stats] = signrank(dbCoh(:,ifreq,ivox),0,'tail','right','alpha',alpha);
            testval(ivox,ifreq) = stats.signedrank;
        end
    end

    onoff = hs;
    clear clu total
    [clu, total] = findcluster(onoff,...
        match_conn, match_conn, minnbchan);
    shuf_clu(iit,:,:) = fp_order_clusters(clu,total);
    shuf_total(iit) = total;
    shuf_testval(iit,:,:) = testval;
    shuf_p(iit,:,:) = ps;
    
    toc
end

%%

if true_total>0 %when at least one true cluster exists    
    
    for iclus = 1:true_total
        
        clear a shuf_clu_val
        a = zeros(size(shuf_testval));
%         a(shuf_clu==iclus) = shuf_testval(shuf_clu==iclus); %compare first clu with first, second with second etc
        a(shuf_clu==1) = shuf_testval(shuf_clu==1); %take always the biggest cluster 
        shuf_clu_val = squeeze(sum(sum(a,2),3));
        
        clear trueCoh temp
        true_clu_val = sum(sum(true_testval(true_clu==iclus))); %scalar
        p(iclus) = sum(shuf_clu_val>true_clu_val)/numel(shuf_clu_val);
        if p(iclus) > 0.01
            break
        end
    end
    
elseif sum(shuf_total)== 0  %when no cluster was found it any iteration
    p= nan;
    
else %when only in shuffled conditions clusters were found
    p = 1;
end


outname = sprintf('%sp_cluster_g_freq_%s_test2D',DIROUT,abs_imag);
save(outname,'p','true_clu','true_p','true_total','-v7.3')