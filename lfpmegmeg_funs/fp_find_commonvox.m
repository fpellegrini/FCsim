function [cp, voxID] = fp_find_commonvox
%finds voxels that exist in all subjects 
%index according to sym_pos, not according to original mni_pos! 

% cd ~/Dropbox/Franziska/Data_MEG_Project/

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

clear mni_pos sym_pos
mni_pos = fp_getMNIpos(patientID{1});
[sym_pos, ~] = fp_symmetric_vol(mni_pos);
cp = sym_pos; %commonvox_pos

for id = 2:numel(patientID)
   
    clear mni_pos sym_pos
    mni_pos = fp_getMNIpos(patientID{id});    
    [sym_pos, ~] = fp_symmetric_vol(mni_pos);
    clear b
    b = sym_pos;
  
    o=1;
    for i= 1:size(cp,1)
        %compare current common voxels with voxels from this id 

        if any(b(:,1)==cp(i,1) & b(:,2)==cp(i,2) & b(:,3)==cp(i,3))

            ind(o,:) = cp(i,:);
            o=o+1;
        end
    end 
    
    clear cp
    cp = ind;
    clear ind
end

for id = 1:numel(patientID)
    
    clear mni_pos sym_pos
    mni_pos = round(fp_getMNIpos(patientID{id}));  
    [sym_pos, ~] = fp_symmetric_vol(mni_pos);
    
    for i = 1:size(cp,1)
        voxID{id}(i) = find(sym_pos(:,1)==cp(i,1) & sym_pos(:,2)==cp(i,2) & sym_pos(:,3)==cp(i,3));
    end
    
    
end 


