function fp_surrogate_coh(patientNumber, njack)

cd ~/Dropbox/MEG_Project/Data

if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

if ~exist('njack','var')
    njack = 3;
end

for id = 1:numel(patientID)

    [~,~,~, coh(1,:,:,:)] = fp_timesensor2sourcecoh(patientID{id}, 0);

    for ijack = 1:njack

        [~,~,~,coh(ijack+1,:,:,:)] = fp_timesensor2sourcecoh(patientID{id}, 1);
    end 
    
    COH{id} = coh;
    clear coh 
    
end 

