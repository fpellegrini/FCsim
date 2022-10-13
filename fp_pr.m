function [pm] = fp_pr(cc,iroi_seed,iroi_tar,half_m)

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

%%
if norm(iroi_seed(:)- iroi_tar(:),'fro')==0
    pm = 0;
else
    
    nints = numel(iroi_seed); %number of interactions 
    nroi = size(cc,1); %number of regions 
    
    if half_m %if only the upper triangle of the roi-roi matrix should be considered
        ind_ = find(triu(ones(nroi),1));
    else
        ind_ = find(ones(nroi));
    end
    
    %get ground-truth interaction pattern 
    label= zeros(nroi,nroi);
    for iint = 1:numel(iroi_seed)
        label(iroi_seed(iint),iroi_tar(iint)) = 1;
        if half_m %when only positive DIFFGC is submitted as cc
            label(iroi_tar(iint),iroi_seed(iint)) = 1;
        end
    end
    lab = label(:);
    tr = find(lab(ind_)==1); %index of true interactions 
    
    %sort estimated FC scores and find the rank of the ground-truth
    %interactions 
    [~, idx] = sort(cc(ind_),'descend');   
    for it = 1:numel(tr)%loop over interactions 
        r1(it) = find(idx==tr(it)); %rank 
    end
    
    %% percentile rank 
    
    pm = mean(1 - (r1 ./ numel(lab(ind_)))); %mean over interactions 
    
    %calculate perfect skill and no skill PR 
    for it = 1:nints
        perfectPm(it) = 1 - (it / numel(lab(ind_))); %perfect skill
        noSkillPm(it) = 1 - (numel(lab(ind_))-it+1)/numel(lab(ind_)); %no skill
    end    
    perfectPm = mean(perfectPm);
    noSkillPm = mean(noSkillPm);
    
    %normalize PR with perfect and no skill PRs 
    pm = (pm-noSkillPm)/(perfectPm-noSkillPm);
    
end