
%% forth case 
clear all

load('./processed_bs/bs_results.mat')
smooth_cortex = 0.35;

% signal generation

rng(1)
% number of ROIs in the Desikan-Kiliany Atlas
nroi = length(cortex.Atlas(3).Scouts);

n_trials = 200;
Lepo = 100;
N = n_trials*Lepo;
lag = 5;
fres = 40;
iroi_seed = 11;
iroi_tar = 65;
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
nroi = length(sub_ind_cortex);
nvox = length(ind_cortex); 

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

% p = randn(ni,1);
p=[0.5; 0.5];
p = p/norm(p);

for in = 1:nroi
    L1 = squeeze(L_forward(:,in,:));
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

clear L1 L2 L_forward L_mix L3

%% megmeg pipeline start
%parameters

id_meg_chan = 1:size(signal_sensor,1);
nmeg = numel(id_meg_chan);

id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(signal_sensor,fres,'MT',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);

%leadfield backward model 
L3 = L_save(:, ind_cortex, :);
for is=1:nvox
    clear L2
    L2 = L3(:,is,:);
    
    %remove radial orientation
    clear u s
    [u, s, v] = svd(squeeze(L2));
    L_backward(:,is,:) = u(:,:)*s(:,1:2);
end

%construct source filter
if strcmp(filtertype,'e')
    A = squeeze(mkfilt_eloreta_v2(L_backward));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(filtertype,'d')

    A=zeros(nmeg,ni,nvox,nfreq);
    
    for ifrq = 1:nfreq
        cCS = CS(:,:,ifrq);
        lambda = mean(diag(real(cCS)))/100;
        
        CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
        
        for is=1:nvox %iterate across nodes
            Lloc=squeeze(L_backward(:,is,:));
            A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
        end
    end
    fqA = 1:nfreq; %This filter is frequency specific.
    nfqA = nfreq;
    
    
elseif strcmp(filtertype,'l')
    
    
end



%%
npcs = 2;

for aroi = 1:nroi 

    %project to source level
    clear A_ CSv
    A_ = A(:, :,ind_roi_cortex{aroi},:);
    nvoxroi = size(A_,3);
    A2 = reshape(A_, [nmeg, ni*nvoxroi, nfqA]);
    
    
    for ifq = 1: nfreq
        CSv(:,:,ifq) = squeeze(A2(:,:,fqA(ifq)))' * CS(:,:,ifq)...
            * squeeze(A2(:,:,fqA(ifq)));
%         n = size(csv,2);
%         csv(1:n+1:end) = real(diag(squeeze(csv)));
%         CSv(ifq,:,:)=csv;
    end
    
    %zscoring
    clear ZS CSz
    ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 3))))./diag(squeeze(sum(real(CSv), 3)))));
    for ifreq = 1:nfreq
        CSz(ifreq,:, :) = ZS'*squeeze(CSv(:, :,ifreq))*ZS;
    end
    
    clear CSs v v5 in V_ D_
    CSs = squeeze(sum(CSz,1)); %covariance
    [V_, D_] = eig(real(CSs));
    [D_, in] = sort(real(diag(D_)), 'descend');
    
    V{aroi} = V_(:,in(1:npcs)); %nregionvoxels*2 x npcs
        
    %     %concatenate filters
    for ifq = 1:nfqA
        P(:, :, aroi,ifq) = A2(:,:,fqA(ifq)) * ZS * real(V{aroi});
    end
end


%apply all filters
CSroi = [];
for ifreq = 1:nfreq
    CSroi(:, :, ifreq) = reshape(P(:,:,:,fqA(ifreq)), nmeg, [])'*CS(:, :, ifreq)...
        *reshape(P(:,:,:,fqA(ifreq)), nmeg, []);
end

%divide by power to obtain coherence
clear Cohroi
for ifreq = 1: fres
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end
   

%%
a = [];

ic=1;
for iroi = 1:nroi 
    
    
    jc=1;
    for jroi = 1:nroi
        
        for ifq = 1:nfqA
            cs_red=[];
            cs_red{1} = Cohroi(ic:ic+npcs-1,ic:ic+npcs-1,ifq);
            cs_red{2} = Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,ifq);
            cs_red{3} = Cohroi(jc:jc+npcs-1,jc:jc+npcs-1,ifq);
            
            caainv=inv(real(cs_red{1})+regu*eye(npcs)*mean(diag(real(cs_red{1}))));
            cab=imag(cs_red{2});
            cbbinv=inv(real(cs_red{3})+regu*eye(npcs)*mean(diag(real(cs_red{3}))));
            X=cab*cbbinv*cab';
            % MIM Ewald Eq. 14
            mim1(iroi,jroi,ifq)=(trace(caainv*X));
            caainvsqrt=sqrtm(caainv);
            Y=caainvsqrt*X*caainvsqrt; %Eq. 23
            [~,s,~]=svd(Y);
            % MIC
            mic1(iroi,jroi,ifq)=sqrt(s(1,1));
        end
        jc = jc+npcs;
    end
    ic=ic+npcs;
end

%%
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

%%
a1 = zeros(size(cortex.Vertices,1),1); 
for ir = 1:nroi 
    a1(ind_roi{ir}) = mm(ir);
end 
load cm17
pos = cortex.Vertices;

xx = zeros(size(a1));
xx([ind_roi{iroi_seed}; ind_roi{iroi_tar}])=0.2;

data_in=xx;
allplots_cortex_BS(cortex, data_in, [min(data_in) max(data_in)],...
    cm17a,'.', smooth_cortex,['test']);
clear data_in

data_in = a1;
allplots_cortex_BS(cortex, data_in, [min(data_in) max(data_in)],...
    cm17a,'.', smooth_cortex,['test1']);
