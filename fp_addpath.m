function fp_addpath

restoredefaultpath
rehash toolboxcache

addpath('~/Dropbox/Master/Masterarbeit/MasterThesis/bb_code/')
addpath('~/Dropbox/Master/Masterarbeit/MasterThesis/wjn_toolbox/')
addpath('~/Dropbox/Master/Masterarbeit/MasterThesis/wjn_toolbox/spm12/')
addpath(genpath('~/Dropbox/Master/Masterarbeit/MasterThesis/wjn_toolbox/spm12/aal'))
addpath(genpath('~/Dropbox/Master/Masterarbeit/MasterThesis/wjn_toolbox/spm12/DAiSS-master'))
spm eeg
close all
addpath(genpath('~/Dropbox/Master/Masterarbeit/MasterThesis/data/'))
addpath(genpath('~/Dropbox/Master/Masterarbeit/MasterThesis/figures/'))
addpath(genpath('~/Dropbox/MEG_Project'))
addpath(genpath('~/Dropbox/Data_MEG_Project'))
