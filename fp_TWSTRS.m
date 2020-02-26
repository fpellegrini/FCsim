function fp_TWSTRS

filtertype = 'd';
abs_imag = 'imag';

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[~, voxID] = fp_find_commonvox;

scores = db_twstrs();
alpha = 0.05;


nsubs = numel(patientID);
nchunk = 50;
%%
o = 1; 
for id = 1:nsubs    
    
    if ~isnan(scores{1,id}{1,2})
        fprintf('Working on subject %s \n',patientID{id})
        
        t_score(o) = scores{1,id}{1,2}; 

        %get neighbouring nodes and node positions
        clear flip_id noEq

        %get flip id and symmetric head
        [flip_id, noEq] = fp_get_flip_id(patientID{id},voxID{id});

        for ichunk = 1:nchunk
            fprintf('Working on chunk %d \n',ichunk)

            if strcmp(filtertype, 'e')
                inname = sprintf('Coherences_e_Patient%s_chunk%d.mat',patientID{id},ichunk); 
            elseif strcmp(filtertype, 'e2D')
                inname = sprintf('Coherences_e2D_Patient%s_chunk%d.mat',patientID{id},ichunk);
            elseif strcmp(filtertype,'d')
                inname = sprintf('Coherences_Patient%s_chunk%d.mat',patientID{id},ichunk); 
            else 
                error('Wrong filtertype!')
            end
            COH1(o,ichunk,:,:,:) = fp_get_coh(inname, noEq, voxID{id},flip_id, abs_imag);

        end 
        o = o+1;
    end
end

[nsubs, nchunk,niit,nfreq,ns] = size(COH1);
COH = reshape(COH1,[nsubs,nchunk*niit,nfreq,ns]);
nit = size(COH,2)-1; %update nit to number of *shuffled* it

%true and shuffled coherences
tCoh = squeeze(COH(:,1,:,:));
sCoh = COH(:,2:end,:,:);

%% true

for ifreq = 1:nfreq 
    for ivox = 1:ns 
        [true_rho(ivox,ifreq), true_pval(ivox, ifreq)] = corr(t_score', tCoh(:,ifreq,ivox), 'Type','Spearman');
    end 
end 

onoff = true_pval < alpha;
kron_conn = fp_get_kron_conn(patientID{id}, voxID{id}, nfreq); %same for every id

%pos
onoff1 = onoff;
onoff1(true_rho<0)=0;
[true_clu_pos, true_total_pos] = fp_get_cluster_components(onoff1,kron_conn);

%neg
onoff1 = onoff;
onoff1(true_rho>0)=0;
[true_clu_neg, true_total_neg] = fp_get_cluster_components(onoff1,kron_conn);

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
