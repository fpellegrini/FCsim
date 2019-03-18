function neighb = fp_find_neighbours(patientNumber)

cd ~/Dropbox/MEG_Project/Data

load(sprintf('BF_Patient%s.mat',patientNumber));
c=sources.grid.pos;

if mean(abs(c(:,1)),1)<1
    c = c.*1000; %avoid numerical errors 
end

dist = pdist(c); %eucledean distance
dist1 = squareform(dist);
dist1(dist1==0)=inf; %set distance at node itself to inf
a=min(dist1);
gridsiz=a(1);

neighb = zeros(size(dist1));

for inode =1:size(neighb,1)
    o = dist1(inode,:);
    
    neighb(inode,find(o<ceil(gridsiz)+2))=1;
    clear o
    
end

% %check wether only outer gridpoints have fewer than 6 neighbours
% outerpoints = find(sum(neighb,1)<6); 
% scatter3(c(:,1),c(:,2),c(:,3))
% hold on
% scatter3(c(outerpoints,1),c(outerpoints,2),c(outerpoints,3),'r+')


