function leadfield = fp_leadfield(patientNumber, changeFilePaths, calculateBF, freqBand, refChannel, plotOn,inode, idir)

%inputs: 
%patientNumber (string), changeFilePaths (0 or 1), calculateBF (0 or 1), freqBand ([low high]),
%refChannel (string with channel name), plotOn (0 or 1),inode (node
%number), idir (1 2 or 3)

cd ~/Dropbox/Data_MEG_Project/
% DIROUT = '~/Dropbox/Data_MEG_Project/figures/topo_leadfields/';
% if ~exist(DIROUT); mkdir(DIROUT); end

if strcmp(patientNumber,'all')
    %change filenames
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'18';'20';'22';'25'}; %'12' has too few sensors
else
    patientID{1} = patientNumber; 
end 
if nargin < 2
    changeFilePaths= 0; %0 or 1
end 
if nargin < 3
    calculateBF = 1; %0 or 1
end
if nargin < 4
    freqBand = [13 30];
end 
if nargin < 5
    refChannel = 'LFP_R01';
end
if nargin < 6
    plotOn = 0;
end 
if nargin < 7
    inode = 1917; 
end 
if nargin< 8
    idir = 1;
end

for id = 1:numel(patientID)
    
    fileName = sprintf('redPLFP%s_off', patientID{id});
    
    clear D
    if changeFilePaths == 1
        %rename filepaths in D
        D = wjn_meg_correct_mri(fileName);
    else
        D = spm_eeg_load(fileName);
    end

    if calculateBF == 1   
        
        %dics according to Neumann 2015
        %alternatively open DAiSS toolbox --> DICS coherence (delete head model module)
        %and choose filename,freqrange,refchannel manually
        fp_dics_coherence(fileName,freqBand, refChannel)

        load('BF.mat')
%         save(sprintf('%s/BF_Patient%s',pwd,patientID{id}),'sources')
        save(sprintf('%s/BF_Patient%s_1',pwd,patientID{id}),'sources', ...
            'inverse','data','features','output','write','-v7.3')
        !rm ./BF.mat
        
    else
        load(sprintf('%s/BF_Patient%s1',pwd,patientID{id}))
    end
    
    clear leadfield
    leadfield = sources.L.MEG;
     
    if plotOn == 1
        %extract channel positions
        D_ft = ftraw(D);
        chanpos = D_ft.grad.chanpos;
        chanpos(D.badchannels,:)=[];
        
        %select leadfield at specific node and direction
        lf = leadfield{inode}(:,idir);
        
        loc = mk_sensors_plane(chanpos(:,[2 1 3]));
        close all
        figure
        showfield_general(lf,loc);
        colormap('jet')
%         caxis([-(10^-12) 10^-12])
%         title(sprintf('Patient %s, node %d, direction %d', patientID{id}, inode, idir))
    
        outname = sprintf('%slf_patient%s_node%d_dir%d.png',DIROUT,patientID{id},inode,idir);
        print(outname,'-dpng');
        close all
    end 


    clear loc fileName D D_ft sources chanpos lf
end  