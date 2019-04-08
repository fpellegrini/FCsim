function fp_dics_coherence(fileName, freqRange, refChannel)
%input: filename, freqrange, refchannel 
%Example: fp_dics_coherence('redPLFP04_off', [13 30], 'LFP_R01')

matlabbatch{1}.spm.tools.beamforming.data.dir = {'/Users/franziskapellegrini/Dropbox/Data_MEG_Project/'};
matlabbatch{1}.spm.tools.beamforming.data.D = {sprintf('/Users/franziskapellegrini/Dropbox/Data_MEG_Project/%s.mat',fileName)};
matlabbatch{1}.spm.tools.beamforming.data.val = 1;
matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{1}.spm.tools.beamforming.data.space = 'MNI-aligned';
matlabbatch{1}.spm.tools.beamforming.data.overwrite = 0;
matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.grid.resolution = 10;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.grid.space = 'MNI template';
matlabbatch{2}.spm.tools.beamforming.sources.plugin.grid.constrain = 'iskull';
matlabbatch{2}.spm.tools.beamforming.sources.visualise = 0;
matlabbatch{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{3}.spm.tools.beamforming.features.woi = [-Inf Inf];
matlabbatch{3}.spfrem.tools.beamforming.features.modality = {'MEG'};
matlabbatch{3}.spm.tools.beamforming.features.fuse = 'no';
matlabbatch{3}.spm.tools.beamforming.features.plugin.csd.foi = freqRange;
matlabbatch{3}.spm.tools.beamforming.features.plugin.csd.taper = 'dpss';
matlabbatch{3}.spm.tools.beamforming.features.plugin.csd.keepreal = 0;
matlabbatch{3}.spm.tools.beamforming.features.plugin.csd.hanning = 0;
matlabbatch{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 0.01;
matlabbatch{3}.spm.tools.beamforming.features.bootstrap = false;
matlabbatch{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.inverse.plugin.dics.fixedori = 'yes';
matlabbatch{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.reference.refchan.name = refChannel;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.reference.refchan.shuffle = 0;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.powmethod = 'lambda1';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.whatconditions.all = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.sametrials = false;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.woi = [-Inf Inf];
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.contrast = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.logpower = false;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.foi = freqRange;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.taper = 'dpss';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.result = 'singleimage';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.scale = 'no';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.modality = 'MEG';
matlabbatch{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.write.plugin.nifti.normalise = 'separate';
matlabbatch{6}.spm.tools.beamforming.write.plugin.nifti.space = 'mni';

matlabbatch{7}.spm.spatial.smooth.data(1) = cfg_dep('Write: Output files', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{7}.spm.spatial.smooth.fwhm = [12 12 12];
matlabbatch{7}.spm.spatial.smooth.dtype = 0;
matlabbatch{7}.spm.spatial.smooth.im = 1;
matlabbatch{7}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch)
