function clu_sorted = fp_order_clusters(clu,total) 

x = hist(clu(:),0:total);
[~, order] = sort(x(2:end),'descend');

clu_sorted = zeros(size(clu)); 

for iclu = 1:total
    clu_sorted(clu==order(iclu))=iclu; 
end