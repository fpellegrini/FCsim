function fp_thresh_free_test

DIRDAT = '/datasabzi/Franziska/';
fwf_s = [];
filtertype ='d';
abs_imag = 'imag';


j_s = 'allsubs';
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[~, voxID] = fp_find_commonvox;


nsubs = numel(patientID);
nchunk = 50;
minnbchan = 3;

alphs = [0.001, 0.002, 0.003, 0.005, 0.007, 0.009, 0.0125, 0.016, 0.021,...
    0.027, 0.034, 0.042, 0.051];
tcrits = [64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52];
%%

for id = 1:nsubs   
    fprintf('Working on subject %s \n',patientID{id})
    
    %get neighbouring nodes and node positions
    clear flip_id noEq
    
    %get flip id and symmetric head
    [flip_id, noEq] = fp_get_flip_id(patientID{id},voxID{id});
    
    for ichunk = 1:nchunk
        fprintf('Working on chunk %d \n',ichunk)
        
        if strcmp(filtertype, 'e')
            inname = sprintf('%sCoherences_e_Patient%s_chunk%d.mat',DIRDAT,patientID{id},ichunk); 
        elseif strcmp(filtertype, 'e2D')
            inname = sprintf('%sCoherences_e2D_Patient%s_chunk%d.mat',DIRDAT,patientID{id},ichunk);
        elseif strcmp(filtertype,'d')
            inname = sprintf('%sCoherences_Patient%s_chunk%d.mat',DIRDAT,patientID{id},ichunk); 
        else 
            error('Wrong filtertype!')
        end
        COH1(id,ichunk,:,:,:) = fp_get_coh(inname, noEq, voxID{id},flip_id, abs_imag);
        
    end   
end

[nsubs, nchunk,niit,nfreq,ns] = size(COH1);
COH = reshape(COH1,[nsubs,nchunk*niit,nfreq,ns]);
nit = size(COH,2)-1; %update nit to number of *shuffled* it

%true and shuffled coherences
tCoh = squeeze(COH(:,1,:,:));
sCoh = COH(:,2:end,:,:);

%% true 
   
kron_conn = fp_get_kron_conn(patientID{id}, voxID{id}, nfreq); %same for every id
conn = fp_find_neighbours(patientID{id});
match_conn = conn(voxID{id},voxID{id});
blobb = zeros(numel(alphs), nfreq,ns); 

o = 1;
for ia = alphs    
    clear onoff
    if o ==1 % get the test statistics
        [~,onoff,true_val] = fp_get_signrank_results(tCoh,sCoh,ia);    
    else 
       onoff = true_val>= tcrits(o);
    end
    
    [true_clu, ~] = fp_get_cluster_findclusters(onoff, match_conn,minnbchan);   
%     [true_clu, ~] = fp_get_cluster_components(onoff,kron_conn);
    cblobb = zeros(size(true_clu)); 
    for iclus = 1:10
        siz = sum(sum(true_clu==iclus));
        cblobb(true_clu==iclus) = siz*tcrits(o);
    end
    blobb(o,:,:) = cblobb; 
    o = o+1; 
end

blobb1 = squeeze(sum(blobb,1)); 
max_blobb = max(blobb1(:)); 
% [blobb_ind1,blobb_ind2] = find(blobb1 == max_blobb); 

%% shuffled

sblobb = zeros(nit, numel(alphs), nfreq,ns); 

for iit = 1:nit
    
    fprintf('Working on iteration %d \n',iit)
    clear onoff cCoh csCoh
    
    %select current shuffled coh
    cCoh = squeeze(sCoh(:,iit,:,:)); 

        
        if iit == 1
            csCoh = sCoh(:,2:end,:,:);
        elseif iit == nit
            csCoh = sCoh(:,1:end-1,:,:);
        else
            csCoh = sCoh(:,[1:iit-1 iit+1:end],:,:);
        end

        
        o = 1;
        for ia = alphs
            clear onoff
            if o ==1 % get the test statistics
                [~,onoff,shuf_val] = fp_get_signrank_results(cCoh,csCoh,ia);
            else
                onoff = shuf_val>= tcrits(o);
            end
            
            [shuf_clu, ~] = fp_get_cluster_findclusters(onoff, match_conn,minnbchan);
            %             [shuf_clu, ~] = fp_get_cluster_components(onoff,kron_conn);
            
            cblobb = zeros(size(shuf_clu));
            for iclus = 1:10
                siz = sum(sum(shuf_clu==iclus));
                cblobb(shuf_clu==iclus) = siz*tcrits(o);
            end
            sblobb(iit,o,:,:) = cblobb;
            o = o+1;
        end
       
end

sblobb1 = squeeze(sum(sblobb,2));
s_max_blobb = squeeze(max(max(sblobb1,[],2),[],3));

%%

p = sum(s_max_blobb>max_blobb)/nit;


%%
outname = sprintf('./TFCE_%s_3minnbchan',filtertype);

save(outname,'p','true_clu','true_val','blobb1', '-v7.3')
