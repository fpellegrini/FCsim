
clear all

load('./processed_bs/bs_results.mat')
smooth_cortex = 0.35;

%% signal generation

rng(1)
% number of ROIs in the Desikan-Kiliany Atlas
nroi = length(cortex.Atlas(3).Scouts);

n_trials = 200;
Lepo = 100;
N = n_trials*Lepo;
lag = 5;
fres = 40;
iroi_seed = 10;
iroi_tar = 50;
filtertype= 'd';
regu=.000001;

s1 = randn(nroi, N);
s1(iroi_tar, :) = circshift(s1(iroi_seed, :), lag, 2);

%SNR = 0.5
s1([iroi_seed iroi_tar],:) = s1([iroi_seed iroi_tar],:)./ norm(s1([iroi_seed iroi_tar],:),'fro');
s1([(1:iroi_seed-1), (iroi_seed+1):(iroi_tar-1),(iroi_tar+1):end],:) = s1([(1:iroi_seed-1),...
    (iroi_seed+1):(iroi_tar-1),(iroi_tar+1):end],:)./ ...
    norm(s1([(1:iroi_seed-1),(iroi_seed+1):(iroi_tar-1),(iroi_tar+1):end],:),'fro');

signal_source = reshape(s1, nroi, Lepo, n_trials);

% ROI labels
labels = {cortex.Atlas(3).Scouts.Label};
%roi inds 
ind_cortex = [];
sub_ind_cortex = [];
ind_roi = {};
sub_ind_roi = {};
for iROI = 1:nroi
  ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices;
  ind_cortex = cat(1, ind_cortex, ind_roi{iROI});
  [~, ind_roi_cortex{iROI}, ~] = intersect(ind_cortex, ind_roi{iROI}); %index of roi voxels in ind_cortex
  sub_ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices(1); 
  sub_ind_cortex = cat(1,sub_ind_cortex, sub_ind_roi{iROI}); 
  [~,sub_ind_roi_cortex{iROI},~] =  intersect(sub_ind_cortex, sub_ind_roi{iROI});%only one voxel per region 
  
end
nvox = length(sub_ind_cortex);
leadfield = leadfield(:, sub_ind_cortex, :);
L = leadfield;
L1 = L;
clear L
for is=1:nroi
    clear L2
    L2 = L1(:,is,:);
    
    %remove radial orientation
    clear u s
    [u, s, v] = svd(squeeze(L2));
    L(:,is,:) = u(:,:)*s(:,1:2);
end

ni = size(L,3);

% p = randn(ni,1);
p=[0.5; 0.5];
p = p/norm(p);

for in = 1:nroi
    L1 = squeeze(L(:,in,:));
    L_mix(:,in) = L1*p;
end

%project to sensors and add white noise
for itrial = 1:n_trials
    clear sig whitenoise noise
    sig = L_mix * signal_source(:,:,itrial);
    sig = sig ./ norm(sig, 'fro');
    
    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal_sensor(:,:,itrial) = sig ./ norm(sig, 'fro');
end

% clear x1 x2 xx whitenoise sig X L1 noise

%% megmeg pipeline start
%parameters

id_meg_chan = 1:size(signal_sensor,1);
nmeg = numel(id_meg_chan);

id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(signal_sensor,fres,'MT',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);

%construct source filter
if strcmp(filtertype,'e')
    A = squeeze(mkfilt_eloreta_v2(L));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(filtertype,'d')

    A=zeros(nmeg,ni,nroi,nfreq);
    
    for ifrq = 1:nfreq
        cCS = CS(:,:,ifrq);
        lambda = mean(diag(real(cCS)))/100;
        
        CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
        
        for is=1:nroi %iterate across nodes
            Lloc=squeeze(L(:,is,:));
            A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
        end
    end
    fqA = 1:nfreq; %This filter is frequency specific.
    nfqA = nfreq;
    
    
elseif strcmp(filtertype,'l')
    
    
end

A2 = reshape(A,nmeg,ni*nroi,nfreq);

for ifq = 1: nfreq
    CSv(:,:,ifq) = squeeze(A2(:,:,fqA(ifq)))' * CS(:,:,ifq)...
        * squeeze(A2(:,:,fqA(ifq)));
end

%%
% CSs = sum(CSv,3); 
% [v,d,~] = eig(real(CSs));
% [~, in] = sort(diag(d),'descend'); 
% 
% V = v(:, in(1:nroi));

% clear CSroi V
% npcs = 5; 
% V= zeros(ns_org*2, npcs*2);
% CSs1 = squeeze(sum(CSv(1:(inode_tar-1)*2,1:(inode_tar-1)*2,:),3)); %covariance
% CSs2 = squeeze(sum(CSv((inode_tar*2)-1:end,(inode_tar*2)-1:end,:),3)); %covariance
% [v1, d1, ~] = eig(real(CSs1));
% [v2, d2, ~] = eig(real(CSs2));
% [~, in1] = sort(diag(d1), 'descend');
% [~, in2] = sort(diag(d2), 'descend');
% 
% V(1:(inode_tar-1)*2,1:npcs) = v1(:,in1(1:npcs)); 
% V((inode_tar*2)-1:end,npcs+1:end) = v2(:,in2(1:npcs)); 

% CSroi = [];
% for ifreq = 1:fres
%     CSroi(:, :, ifreq) = V'*CSv(:, :, ifreq)*V;
% end

CSroi = CSv;
clear Cohroi

%divide by power to obtain coherence
for ifreq = 1: fres
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end
   

%%
chan = ni; 

for iroi = 1:nroi 
    
    ic = ((iroi-1)*2)+1:((iroi-1)*2)+2;
    for jroi = 1:nroi
        jc = ((jroi-1)*2)+1:((jroi-1)*2)+2;
    
        for ifq = 1:nfqA
            cs_red=[];
            cs_red{1} = Cohroi(ic,ic,ifq); %Caa
            cs_red{2} = Cohroi(ic,jc,ifq); %Cab
            cs_red{3} = Cohroi(jc,jc,ifq); %Cbb

            caainv=inv(real(cs_red{1})+regu*eye(chan)*mean(diag(real(cs_red{1}))));
            cab=imag(cs_red{2});
            cbbinv=inv(real(cs_red{3})+regu*eye(chan)*mean(diag(real(cs_red{3}))));
            X=cab*cbbinv*cab';
            % MIM Ewald Eq. 14
            mim1(iroi,jroi,ifq)=(trace(caainv*X));
            caainvsqrt=sqrtm(caainv);
            Y=caainvsqrt*X*caainvsqrt; %Eq. 23
            [~,s,~]=svd(Y);
            % MIC
            mic1(iroi,jroi,ifq)=sqrt(s(1,1));
        end
    end
end

mic = sum(mic1,3);
mim = sum(mim1,3);

mc = sum(mic,2); 
mm = sum(mim,2); 

imagesc(mic)
figure
imagesc(mim)
figure;
plot((mc- mean(mc))./std(mc(:)))
hold on 
plot((mm - mean(mm))./std(mm(:)))
legend('mic','mim')
grid on 

