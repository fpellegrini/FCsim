function [clu, total] = fp_get_cluster_components_megmeg(onoff,kron_conn)
% keyboard
[nroi,nroi,nfreq] = size(onoff);

u = onoff(:); %should be the same indexing like in kron_conn now; nkron x 1
ind = find(u==1); %remember indeces of super-threshold coherences
A = kron_conn;
A(u==0,:)=[]; %pass the neighbourhood structure only for the super-threshold voxels
A(:,u==0)=[];

%components assigns every voxel to a cluster, even if this means that every voxel is its own cluster
clear ci x clu
[ci, x] = components(A); %x is the histogram of the clusters
clu = zeros(size(kron_conn,1),1);%refill with sub-threshold voxels
clu(ind)= ci; %nkron x 1
clu = reshape(clu,[nroi, nroi, nfreq]); 

total = numel(x);
clu = fp_order_clusters(clu,total);