function [signal_sensor, gt,L_save,iroi_seed,iroi_tar] = fp_generate_mim_signal...
    (params, fres,n_trials, D)

iroi_seed = randi(D.nroi,params.iInt,1);
iroi_tar = randi(D.nroi,params.iInt,1);

%set parameters
Lepo = 100;
N = n_trials*Lepo;
lag = randi(20,params.iInt,params.iReg);
id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;

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
s1(iroi_s,:,:) = s1(iroi_s,:,:)./ norm(reshape(s1(iroi_s,:,:),numel(iroi_s),[]),'fro');
s1(iroi_ns,:,:) = s1(iroi_ns,:,:)./ norm(reshape(s1(iroi_ns,:,:),numel(iroi_ns),[]),'fro');

%generate ground truth imaginary coherence
signal_gt = reshape(permute(s1,[2 1 3]), params.iReg*D.nroi, Lepo, n_trials); %permute that it matches with sub_ind_cortex
CS_gt = fp_tsdata_to_cpsd(signal_gt,fres,'WELCH',...
    1:D.nroi*params.iReg, 1:D.nroi*params.iReg, id_trials_1, id_trials_2);
CS_gt(:,:,1)=[];
for ifreq = 1: fres
    clear pow
    pow = real(diag(CS_gt(:,:,ifreq)));
    imCoh_gt(:,:,ifreq) = abs(imag(CS_gt(:,:,ifreq)./ sqrt(pow*pow')));
end
gt1 = sum(sum(imCoh_gt,2),3);
gt = squeeze(mean(reshape(gt1,2,[]),1));

%leadfield for forward model
L_save = D.leadfield;
L3 = L_save(:, D.sub_ind_cortex, :);
for is=1:D.nroi*params.iReg
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

for in = 1 : D.nroi*params.iReg
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
    sig = params.isnr*sig + (1-params.isnr)*whitenoise;
    signal_sensor(:,:,itrial) = sig ./ norm(sig, 'fro');
end