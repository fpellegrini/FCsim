function neighb = fp_find_neighbours(patientNumber)

mni_pos = fp_getMNIpos(patientNumber);
[c, ~] = fp_symmetric_vol(mni_pos);

dist = pdist(c); %eucledean distance
dist1 = squareform(dist);
dist1(dist1==0)=inf; %set distance at node itself to inf
a=min(dist1);
gridsiz=a(1);

neighb = zeros(size(dist1));

for inode =1:size(neighb,1)
    o = dist1(inode,:);
    
    neighb(inode,find(o==gridsiz))=1;
    clear o
    
    neighb(inode,inode) = 1; %voxel is also neighbour to itself
    
end
% 
% %check visually wether only outer gridpoints have fewer than 6 neighbours
% outerpoints = find(sum(neighb,1)<6); 
% figure
% scatter3(c(:,1),c(:,2),c(:,3))
% hold on
% scatter3(c(outerpoints,1),c(outerpoints,2),c(outerpoints,3),'r+')
% pause(2)


