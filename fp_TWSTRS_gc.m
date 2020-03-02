function fp_TWSTRS_gc

DIROUT = './';
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[~, voxID] = fp_find_commonvox;

scores = db_twstrs();
alpha = 0.05;
%%
load(sprintf('%sDIFFGC.mat',DIROUT));
[nsubs,nvox,nside,nfreq] = size(DIFFGC);
%%
o = 1;
for id = 1:nsubs
    if ~isnan(scores{1,id}{1,2})
        fprintf('Working on subject %s \n',patientID{id})
        
        t_score(o) = scores{1,id}{1,2};
        diffgc_t(o,:,:,:) = DIFFGC(id,:,:,:);
        o = o+1;
        
    end
end
clear o
%% true

for iside = 1:2 
    for ifreq = 1:nfreq 
        for ivox = 1:nvox 
            [true_rho(ivox,ifreq, iside), true_pval(ivox, ifreq, iside)] ...
                = corr(t_score', diffgc_t(:,ivox,iside,ifreq), 'Type','Spearman');
        end 
    end 
end

kron_conn = fp_get_kron_conn(patientID{id}, voxID{id}, nfreq); %same for every id
for iside = 1:2 
    clear onoff onoff1 
    onoff = true_pval(:,:,iside) < alpha;

    %pos
    onoff1 = onoff;
    onoff1(true_rho(:,:,iside)<0)=0;
    [true_clu_pos(:,:,iside), true_total_pos(:,iside)] = fp_get_cluster_components(onoff1,kron_conn);

    %neg
    onoff1 = onoff;
    onoff1(true_rho(:,:,iside)>0)=0;
    [true_clu_neg(:,:,iside), true_total_neg(:,:,iside)] = fp_get_cluster_components(onoff1,kron_conn);
end

%% shuffled

for iit = 1:nit
    fprintf('Working on iteration %d \n',iit);
    clear pval_shuf onoff onoff1 
    
    for ifreq = 1:nfreq
        for ivox = 1:ns           
            [rho_shuf(iit,ivox,ifreq), pval_shuf(ivox, ifreq)] = corr(t_score', sCoh(:,iit,ifreq,ivox), 'Type','Spearman');
        end
    end    
    onoff = pval_shuf < alpha;
    
    %pos 
    onoff1 = onoff;
    onoff1(rho_shuf(iit,:,:,:)<0)=0;
    [shuf_clu(iit,:,:,1), shuf_total(iit,1)] = fp_get_cluster_components(onoff1,kron_conn);
    
    %neg
    onoff1 = onoff;
    onoff1(rho_shuf(iit,:,:,:)>0)=0;
    [shuf_clu(iit,:,:,2), shuf_total(iit,2)] = fp_get_cluster_components(onoff1,kron_conn);
end



%%

p_pos = fp_get_cluster_p_gc_new(true_total_pos, shuf_total, true_rho, rho_shuf, true_clu_pos, shuf_clu, 0);
p_neg = fp_get_cluster_p_gc_new(true_total_neg, shuf_total, true_rho, rho_shuf, true_clu_neg, shuf_clu, 0);

outname = './TWSTRS.mat';
save(outname,'p_pos','p_neg','true_total_pos','true_total_neg','true_clu_pos','true_clu_neg','true_pval','true_rho','-v7.3')
