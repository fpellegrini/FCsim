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
    load(sprintf('BF_Patient%s.mat',patientID{id}))

    a = sources.grid.pos(:,1);
%     hist(a,30)
%     xlabel('position in x direction') 
%     ylabel('number of grid points')
%     
    b = -a;
%     plot(a,b,'.')
%     grid on

    dist = pdist(c);
    dist_ = squareform(dist);
    imagesc(dist_)


    
    c=sources.grid.pos;
    d=c;
    d(:,1)=-c(:,1);
    
    scatter3(c(:,1),c(:,2),c(:,3),'.')
    hold on 
    scatter3(d(:,1),d(:,2),d(:,3),'r.')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    
    clear outname
    outname = sprintf('%shist_Patient%s.png',DIROUT,patientID{id});
    print(outname,'-dpng');
    close all
end 
