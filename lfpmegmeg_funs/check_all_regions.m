a = roi_id;
a(a==0)=[];
a1= sort(a);

count = zeros(1,nroi);
j = 1;
count(j) = 1; 

for ii = 2:numel(a1) 
    if a1(ii)~=a1(ii-1)
        j= j+1;
    end
    count(j) = count(j)+1;
end
    
bad = find((count<2) & (count~=0));
rois = u_roi_id; 
rois(1) = [];
bad_id = rois(bad);

load('ROI_MNI_V5_List.mat')

for ii = 1:120 
    c(ii) = ROI(ii).ID;
end

for ii=1:numel(bad_id)
    f = find(c==bad_id(ii));
    bad_labels{ii,:} = ROI(f).Nom_L;
end