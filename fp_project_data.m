function sData = fp_project_data(patientNumber)

cd ~/Dropbox/MEG_Project/Data

if nargin < 1
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end  

for id = 1: numel(patientID)
    
    load(sprintf('BF_Patient%s.mat',patientID{id}));

    nNodes = numel(inverse.MEG.W);
    nChans = numel(inverse.MEG.W{1});
    nTrials = data.D.ntrials;
    nSamples = data.D.nsamples; %samples per trial

    %filters
    filter = nan(nNodes,nChans);
    for iNode = 1: nNodes
        filter(iNode,:) = inverse.MEG.W{iNode};
    end 

    %data
    dat = data.D(:,:,:);
    dat(126:end, :,:) = [];%select only MEG channels
    dat(data.D.badchannels,:,:)=[];

    %data on source level
    sData = nan(nNodes,nSamples,nTrials);

    %project data through filter
    for iTrial = 1: nTrials
        for iSample = 1: nSamples 
            sData(:,iSample,iTrial) = dat(:,iSample,iTrial)' * filter';
        end 
    end 
    
    outname = sprintf('%s/sData_Patient%s',pwd,patientID{id});
    save(outname,'sData','-v7.3')
    clearvars -except patientID id
end


