function [sig,brain_noise,sensor_noise, gt,L_save,iroi_seed,iroi_tar,D, fres, n_trials] = fp_generate_mim_signal...
    (params,D,DIROUT1)

%if second condition of lag, then load parameters of first condition 
flag = true; 
if params.ip==6 
    params_save = params; 
    load(sprintf('%smim_lag/%d.mat',DIROUT1,params.iit));
    params = params_save; 
    clear params_save
    flag = false; 
end

if flag
    iroi_seed = randi(D.nroi,params.iInt,1);
    iroi_tar = randi(D.nroi,params.iInt,1);
    while iroi_seed==iroi_tar
        iroi_tar = randi(D.nroi,params.iInt,1);
    end
end

%set parameters
Lepo = 50;
fres = Lepo/2;
n_trials = 400;
N = n_trials*Lepo;
id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;

%set random small or large lag
if params.ilag == 1
    lag = randi([0, 5],params.iInt,params.iReg);
else
    lag = randi([5, 20],params.iInt,params.iReg);
end


if flag
    %random signal generation and shifting 
    s1 = randn(D.nroi,params.iReg, N);
end

%save this state of s1 for later use 
if params.ip==1
    s1_save = s1;
end

for iint = 1:params.iInt
    for ireg = 1:params.iReg
        s1(iroi_tar(iint),ireg, :) = circshift(s1(iroi_seed(iint),ireg, :), lag(iint,ireg), 3);
    end
end
%with real data, we would need a normalization step here 
signal_gt = reshape(permute(s1,[2 1 3]), params.iReg*D.nroi, Lepo*n_trials);

%% generate ground truth imaginary coherence
signal_gt_trials = reshape(permute(s1,[2 1 3]), params.iReg*D.nroi, Lepo, n_trials);
CS_gt = fp_tsdata_to_cpsd(signal_gt_trials,fres,'WELCH',...
    1:D.nroi*params.iReg, 1:D.nroi*params.iReg, id_trials_1, id_trials_2);
CS_gt(:,:,[1])=[];
for ifreq = 1: fres
    clear pow
    pow = real(diag(CS_gt(:,:,ifreq)));
    gt1(:,:,ifreq) = CS_gt(:,:,ifreq)./ sqrt(pow*pow');
end
if params.iReg~=1 %mim across two voxels of one region 
    [gt_mic,gt_mim]= fp_mim(gt1,repmat(params.iReg,D.nroi,1));
    gt.mic = gt_mic; 
    gt.mim = gt_mim; 
else 
    gt.mic = abs(imag(gt1));
    gt.mim = abs(imag(gt1)); 
end
clear gt1 gt_mic gt_mim

%% leadfield for forward model

L_save = D.leadfield;
L3 = L_save(:, D.sub_ind_cortex, :);
normals = D.normals(D.sub_ind_cortex,:)'; 
for is = 1:numel(D.sub_ind_cortex)
    L_mix(:,is) = squeeze(L3(:,is,:))*squeeze(normals(:,is));
end

%select signal L and noise L 
sig_ind = [];
for ii = 1:params.iReg
    sig_ind = [sig_ind, (iroi_seed.*params.iReg)-(ii-1), (iroi_tar.*params.iReg)-(ii-1)];
end
L_sig = L_mix(:,sig_ind);
noise_ind = 1:D.nroi*params.iReg;
noise_ind(sig_ind)=[];
L_noise = L_mix(:,noise_ind);

%project to sensors and add white noise 
%signal
sig = L_mix(:,sig_ind) * signal_gt(sig_ind,:);
sig = sig ./ norm(sig, 'fro');    
%brain noise
brain_noise = L_mix(:,noise_ind) * signal_gt(noise_ind,:);
brain_noise = brain_noise ./ norm(brain_noise, 'fro');    
%white noise
if flag
    sensor_noise = randn(size(sig));
    sensor_noise = sensor_noise ./ norm(sensor_noise, 'fro');
end

if params.ip==1
    fprintf('Saving lag stuff... \n')
    dir1 =  sprintf('%smim_lag/',DIROUT1);
    if ~exist(dir1); mkdir(dir1); end
    outname = sprintf('%smim_lag/%d.mat',DIROUT1,params.iit);
    s1 = s1_save;
    save(outname,'iroi_seed','iroi_tar','D','sensor_noise','s1','-v7.3')   
end
