function [sig,brain_noise,sensor_noise, gt,L_save,iroi_seed,iroi_tar] = fp_generate_mim_signal...
    (params, fres,n_trials, D)

iroi_seed = randi(D.nroi,params.iInt,1);
iroi_tar = randi(D.nroi,params.iInt,1);

%set parameters
Lepo = 100;
N = n_trials*Lepo;
id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;

%set random small or large lag

if params.ilag == 1
    lag = randi([0, 5],params.iInt,params.iReg);
else
    lag = randi([5, 20],params.iInt,params.iReg);
end


%random signal generation and shifting 
s1 = randn(D.nroi,params.iReg, N);
for iint = 1:params.iInt
    for ireg = 1:params.iReg
        s1(iroi_tar(iint),ireg, :) = circshift(s1(iroi_seed(iint),ireg, :), lag(iint,ireg), 3);
    end
end

%normalize signal strength
iroi_s = sort(unique([iroi_seed iroi_tar]),'ascend');
iroi_ns = 1:D.nroi; 
iroi_ns(iroi_s)=[];
s1(iroi_s,:,:) = (s1(iroi_s,:,:)./ norm(reshape(s1(iroi_s,:,:),numel(iroi_s),[]),'fro'));
s1(iroi_ns,:,:) = (s1(iroi_ns,:,:)./ norm(reshape(s1(iroi_ns,:,:),numel(iroi_ns),[]),'fro'));

%generate ground truth imaginary coherence
signal_gt = reshape(permute(s1,[2 1 3]), params.iReg*D.nroi, Lepo, n_trials); %permute that it matches with sub_ind_cortex
CS_gt = fp_tsdata_to_cpsd(signal_gt,fres,'WELCH',...
    1:D.nroi*params.iReg, 1:D.nroi*params.iReg, id_trials_1, id_trials_2);
CS_gt(:,:,1)=[];
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
    gt.mic = gt1;
    gt.mim = gt1; 
end
clear gt1 gt_mic gt_mim

%leadfield for forward model
L_save = D.leadfield;
L3 = L_save(:, D.sub_ind_cortex, :);
normals = D.normals(D.sub_ind_cortex,:)'; 
for is = 1:numel(D.sub_ind_cortex)
    L_mix(:,is) = squeeze(L3(:,is,:))*squeeze(normals(:,is));
end

sig_ind = [];
for ii = 1:params.iReg
    sig_ind = [sig_ind, (iroi_seed.*params.iReg)-(ii-1), (iroi_tar.*params.iReg)-(ii-1)];
end
L_sig = L_mix(:,sig_ind);
noise_ind = 1:D.nroi*params.iReg;
noise_ind(sig_ind)=[];
L_noise = L_mix(:,noise_ind);

%project to sensors and add white noise
for itrial = 1:n_trials
    clear sig sensor_noise noise brain_noise noise 
    
    %signal
    sig{itrial} = L_mix(:,sig_ind) * signal_gt(sig_ind,:,itrial);
    sig{itrial} = sig{itrial} ./ norm(sig{itrial}, 'fro');    
    %brain noise
    brain_noise{itrial} = L_mix(:,noise_ind) * signal_gt(noise_ind,:,itrial);
    brain_noise{itrial} = brain_noise{itrial} ./ norm(brain_noise{itrial}, 'fro');    
    %white noise
    sensor_noise{itrial} = randn(size(sig{itrial}));
    sensor_noise{itrial} = sensor_noise{itrial} ./ norm(sensor_noise{itrial}, 'fro');

end