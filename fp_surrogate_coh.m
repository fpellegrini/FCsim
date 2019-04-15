function fp_surrogate_coh(patientNumber, nit, DIROUT, DIRLOG)

if ~exist('DIROUT','var')
    DIROUT =  '~/Dropbox/Data_MEG_Project/';
end
if ~exist(DIROUT); mkdir(DIROUT); end

if ~exist(DIRLOG); mkdir(DIRLOG); end


if isempty(patientNumber)
    patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
else
    patientID{1} = patientNumber;
end

if ~exist('nit','var')
    nit = 1000;
end

for id = 1:numel(patientID)
    fprintf('Working on subject %s',patientID{id})
    logname = sprintf('%s',patientID{id});
    
    if ~exist(sprintf('%s%s_work',DIRLOG,logname)) & ~exist(sprintf('%s%s_done',DIRLOG,logname))
        eval(sprintf('!touch %s%s_work',DIRLOG,logname))
    
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
        eval(sprintf('!mv %s%s_work %s%s_done',DIRLOG,logname,DIRLOG,logname))
    end  
end 

