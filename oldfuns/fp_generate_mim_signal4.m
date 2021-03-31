function [signal_sensor, gt,L_save] = fp_generate_mim_signal4(iroi_seed1,iroi_tar1,iroi_seed2, iroi_tar2, fres, n_trials,sub_ind_cortex,leadfield)

fprintf(['seed1 ' num2str(iroi_seed1) ', tar1 ' num2str(iroi_tar1) '; seed2 ' num2str(iroi_seed2) ', tar2 ' num2str(iroi_tar2) '\n'])
tic

%set parameters
% rng(1)
Lepo = 100;
N = n_trials*Lepo;
lag = 5;
id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
nroi = length(sub_ind_cortex);

%signal generation
s1 = randn(nroi, N);
s1(iroi_tar1, :) = circshift(s1(iroi_seed1, :), lag, 2);
s1(iroi_tar2, :) = circshift(s1(iroi_seed2, :), lag, 2);

iroi_s = sort([iroi_seed1,iroi_tar1, iroi_seed2, iroi_tar2],'ascend');
iroi_ns = 1:68; 
iroi_ns(iroi_s)=[];

%SNR = 0.5
s1(iroi_s,:) = s1(iroi_s,:)./ norm(s1(iroi_s,:),'fro');
s1(iroi_ns,:) = s1(iroi_ns,:)./ norm(s1(iroi_ns,:),'fro');

%generate ground truth imaginary coherence
signal_gt = reshape(s1, nroi, Lepo, n_trials);
CS_gt = fp_tsdata_to_cpsd(signal_gt,fres,'WELCH',1:nroi, 1:nroi, id_trials_1, id_trials_2);
CS_gt(:,:,[1 47:end])=[];
for ifreq = 1: fres
    clear pow
    pow = real(diag(CS_gt(:,:,ifreq)));
    imCoh_gt(:,:,ifreq) = abs(imag(CS_gt(:,:,ifreq)./ sqrt(pow*pow')));
end
gt = sum(sum(imCoh_gt,2),3);


nroi = length(sub_ind_cortex);

L_save = leadfield;

%leadfield for forward model
L3 = L_save(:, sub_ind_cortex, :);
for is=1:nroi
    clear L2
    L2 = L3(:,is,:);
    
    %remove radial orientation
    clear u s
    [u, s, v] = svd(squeeze(L2));
    L_forward(:,is,:) = u(:,:)*s(:,1:2);
end

ni = size(L_forward,3);

p = randn(ni,1);
p = p/norm(p);

for in = 1:nroi
    L1 = squeeze(L_forward(:,in,:));
    L_mix(:,in) = L1*p;
end

%project to sensors and add white noise
for itrial = 1:n_trials
    clear sig whitenoise noise
    sig = L_mix * signal_gt(:,:,itrial);
    sig = sig ./ norm(sig, 'fro');
    
    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal_sensor(:,:,itrial) = sig ./ norm(sig, 'fro');
end