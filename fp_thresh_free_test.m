function fp_thresh_free_test

fwf_s = [];
filtertype ='e2D';
abs_imag = 'imag';


j_s = 'allsubs';
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
[~, voxID] = fp_find_commonvox;


nsubs = numel(patientID);
nchunk = 1;
minnbchan = 2;

%from Gibbons, J. D., and S. Chakraborti. Nonparametric Statistical Inference, 2011.
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
            inname = sprintf('Coherences_e_Patient%s_chunk%d.mat',patientID{id},ichunk); 
        elseif strcmp(filtertype, 'e2D')
            inname = sprintf('Coherences_e2D_Patient%s_chunk%d.mat',patientID{id},ichunk);
        elseif strcmp(filtertype,'d')
            inname = sprintf('Coherences_Patient%s_chunk%d.mat',patientID{id},ichunk); 
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

o = 1;
for ia = alphs    
    clear onoff
    if o ==1 % get the test statistics
        [~,onoff,true_val] = fp_get_signrank_results(tCoh,sCoh,ia);    
    else 
       onoff = true_val>= tcrits(o);
    end
    
    [true_clu, ~] = fp_get_cluster_findclusters(onoff, match_conn,minnbchan);    
    for iclus = 1:10
        icc = sum(sum(true_clu==iclus));
        ict = sum(true_val(true_clu==iclus));
        TFCE(iclus,o) = (icc*ict)/(10^5);
    end
    o = o+1; 
end

tTFCE = sum(TFCE,2);
%% shuffled

for iit = 2:nit 
    
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
                [~,onoff,shuf_val] = fp_get_signrank_results(tCoh,sCoh,ia);
            else
                onoff = shuf_val>= tcrits(o);
            end

            [shuf_clu, ~] = fp_get_cluster_findclusters(onoff, match_conn,minnbchan);
            for iclus = 1:10
                icc = sum(sum(shuf_clu==iclus));
                ict = sum(shuf_val(shuf_clu==iclus));
                TFCE(iclus,o) = (icc*ict)/(10^5);
            end
            o = o+1;
        end
        
        sTFCE(iit,:) = sum(TFCE,2);
    
end

for iclus = 1:10
    p(iclus) = sum(sTFCE(:,1)>tTFCE(iclus))/nit;
end


