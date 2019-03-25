function fp_check_node_symmetry(patientNumber) 

cd ~/Dropbox/MEG_Project/Data
DIROUT = '~/Dropbox/MEG_Project/Data/figures/check_node_symmetry/';
if ~exist(DIROUT); mkdir(DIROUT); end

if ~exist('patientNumber','var')
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; %'12' has too few sensors
else
    patientID{1} = patientNumber;
end

for id = 1:numel(patientID) 

    mni_pos = fp_getMNIpos(patientID{id});
    symm_pos = fp_symmetric_vol(mni_pos);
    flipped = fp_flip_vol(symm_pos);
    
    scatter3(symm_pos(:,1),symm_pos(:,2),symm_pos(:,3),'o')
    hold on 
    scatter3(flipped(:,1),flipped(:,2),flipped(:,3),'r+')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    
    clear outname
    outname = sprintf('%spos_Patient%s.png',DIROUT,patientID{id});
    print(outname,'-dpng');
    close all
end 
