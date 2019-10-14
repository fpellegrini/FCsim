function fp_surrogate_coh_e2D(patientNumber, DIROUT, DIRLOG)

fp_addpath

if ~exist('DIROUT','var')
    DIROUT =  '~/Dropbox/Franziska/Data_MEG_Project/';
end
if ~exist(DIROUT); mkdir(DIROUT); end

if ~exist(DIRLOG); mkdir(DIRLOG); end


if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
else
    patientID{1} = patientNumber;
end

nit = 20;
nchunks = 50;


for id = 1:numel(patientID)
    for ichunk = 1:nchunks
        
        fprintf('Working on subject %s, chunk %d',patientID{id}, ichunk)
        logname = sprintf('%s_%d',patientID{id},ichunk);
        
        if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
            eval(sprintf('!touch %s%s_work',DIRLOG,logname))
            
            if ichunk ==1
                %true coherence
                [coh(1,:,:,:)] = fp_timesensor2sourcecoh_e(patientID{id}, 0);
                for iit = 2:nit
                    iit
                    %shuffled coherences
                    [coh(iit,:,:,:)] = fp_timesensor2sourcecoh_e(patientID{id}, 1);
                end
            else
                for iit = 1:nit
                    iit
                    %only shuffled coherences
                    [coh(iit,:,:,:)] = fp_timesensor2sourcecoh_e(patientID{id}, 1);
                end
            end
            
            
            outname = sprintf('%sCoherences_e2D_Patient%s_chunk%d',DIROUT, patientID{id},ichunk);
            save(outname,'coh','-v7.3')            
            clear coh
            
            eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
        end
    end
end

