function [signal_sensor, gt,L_save] = fp_generate_mim_signal...
    (params, fres,n_trials, D)

iroi_seed = randi(D.nroi,params.iInt,1);
iroi_tar = randi(D.nroi,params.iInt,1);

%set parameters
Lepo = 100;
N = n_trials*Lepo;
lag = randi(20,params.iInt*params.iReg,1);
id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;

%random signal generation and shifting 
s1 = randn(D.nroi, N);
for iint = 1:params.iInt
    s1(iroi_tar(iint), :) = circshift(s1(iroi_seed(iint), :), lag(iint), 2);
end

iroi_s = sort(unique([iroi_seed iroi_tar]),'ascend');
iroi_ns = 1:D.nroi; 
iroi_ns(iroi_s)=[];

%normalize signal strength
s1(iroi_s,:) = s1(iroi_s,:)./ norm(s1(iroi_s,:),'fro');
s1(iroi_ns,:) = s1(iroi_ns,:)./ norm(s1(iroi_ns,:),'fro');

%generate ground truth imaginary coherence
signal_gt = reshape(s1, D.nroi, Lepo, n_trials);
CS_gt = fp_tsdata_to_cpsd(signal_gt,fres,'WELCH',1:D.nroi, 1:D.nroi, id_trials_1, id_trials_2);
CS_gt(:,:,[1 47:end])=[];
for ifreq = 1: fres
    clear pow
    pow = real(diag(CS_gt(:,:,ifreq)));
    imCoh_gt(:,:,ifreq) = abs(imag(CS_gt(:,:,ifreq)./ sqrt(pow*pow')));
end
gt = sum(sum(imCoh_gt,2),3);


%leadfield for forward model
L_save = D.leadfield;
L3 = L_save(:, D.sub_ind_cortex, :);
for is=1:D.nroi
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