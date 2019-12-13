function p = fp_get_cluster_p_gc(true_total, shuf_total, true_val, shuf_val, true_clu, shuf_clu, fwf)
%fwf = testing method (less conservative): test first cluster against
%first, second with second etc. Otherwise clusters are always compared
%against the largest shuffled cluster. 

if (~exist(fwf,'var'))|(fwf~=1)
    fwf = 0;
end 

if true_total>0 %when at least one true cluster exists    
    
    for iclus = 1:true_total
        
        clear a shuf_clu_val
        a = zeros(size(shuf_val));
        
        if fwf == 1
            a(shuf_clu==iclus) = shuf_val(shuf_clu==iclus); %compare first clu with first, second with second etc
        else
            a(shuf_clu==1) = shuf_val(shuf_clu==1); %take always the biggest cluster
        end
        
        shuf_clu_val = squeeze(sum(sum(sum(a,2),3),4));
        
        true_clu_val = sum(sum(sum(true_val(true_clu==iclus)))); %scalar
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