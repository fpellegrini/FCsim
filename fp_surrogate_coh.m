function fp_surrogate_coh(patientNumber, nit)

cd ~/Dropbox/MEG_Project/Data

DIROUT = '~/Dropbox/MEG_Project/Data/';
if ~exist(DIROUT); mkdir(DIROUT); end

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

if ~exist('nit','var')
    nit = 50;
end

for id = 1:numel(patientID)
    
    %true coherence 
    [coh(1,:,:,:)] = fp_timesensor2sourcecoh(patientID{id}, 0);

    for iit = 1:nit
        
        iit
        %shuffled coherences
        [coh(iit+1,:,:,:)] = fp_timesensor2sourcecoh(patientID{id}, 1);
    end 

    outname = sprintf('%sCoherences_Patient%s',DIROUT, patientID{id});
    save(outname,'coh','-v7.3')
    
    clear coh 
    
end 

