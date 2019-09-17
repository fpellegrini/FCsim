clear all
patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'}; 
id = 5; 
D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
load(sprintf('Filter_Patient%s.mat',patientID{id}));%Filter and whole CS
ns= size(A,2);

z = zeros(1,ns);

z(2100)=0.8;

outname = 'pow4.nii';
fp_data2nii(z,id,[],outname)