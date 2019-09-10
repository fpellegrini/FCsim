function roi_conn = fp_roi_conn

DIROUT = '~/Dropbox/Franziska/Data_MEG_Project/';

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id=6;
% [commonvox_pos, voxID] = fp_find_commonvox;
mni_pos = fp_getMNIpos(patientID{id});
[sym_pos, ~] = fp_symmetric_vol(mni_pos);
% match_pos = sym_pos(voxID{id},:);

conn = fp_find_neighbours(patientID{id}); 
% match_conn = conn(voxID{id},voxID{id});
ns = size(conn,1);

for ii = 1: ns
    [~,~,roi_id(ii)]=fp_get_mni_anatomy(sym_pos(ii,:));
end 

u_roi_id = sort(unique(roi_id));
load('ROI_MNI_V5_List.mat')

nroi =length(ROI);
roi_conn = zeros(nroi,nroi);

u=[];
for o = 1:nroi 
    u = [u ROI(o).ID];
end

n1= find(~ismember(u,u_roi_id)); %hinzufuegen
n2 = find(~ismember(u_roi_id,u)); %wegnehmen

u_roi_id(n2) = [];
t = 1;
for iroi =1:nroi 
    if any(ismember(iroi,n1))
        roivox{iroi} = nan;
    else
        roivox{iroi} = find(roi_id == u_roi_id(t));
        t=t+1;
    end
end



for iroi = 1:nroi
    for jroi = 1:nroi
        
        if any(~isnan(roivox{iroi})) && any(~isnan(roivox{jroi}))
            clear a
            a = conn(roivox{iroi},roivox{jroi});

            if any(a(:)>0)
                roi_conn(iroi,jroi) = 1;
            end
        end
    end
end

roi_conn(n1,:)=[];
roi_conn(:,n1)=[];

outname = sprintf('%sroi_conn',DIROUT);
save(outname,'roi_conn','-v7.3')

