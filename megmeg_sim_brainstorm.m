clear all

%% signal generation
%parameters

load('./processed_bs/bs_results.mat')
smooth_cortex = 0.35;

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2;
% nfreq = 46;
delay = 30; %in samples, is equal to 100 ms
inode1 = 56; %randi(size(A,2),1);
inode2 = 949;
n_rand_nodes = 100;
rand_nodes = randi(950,1,n_rand_nodes); %generate 100 nodes with random signals
isens = randi(125,1);
idir = randi(2,1);
alph = 0.5;

%load MEG time series and CS
D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
X = D(:,:,:);
id_meg_chan = 1:125;
X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;
 
%correlated nodes 
x1 = squeeze(X(isens,:,:)); %random time meg series
x2 = x1(delay:end,:); %second time series with time delay
x1 = x1(1:end-delay+1,:);
xx(1,:,:)= x1;
xx(2,:,:) = x2;

%uncorrelated nodes
for ii =1:n_rand_nodes
    ind = randperm(70);
    xx_noise(ii,:,:) = squeeze(X(isens,delay:end,ind));
end

n_trials = size(x1,2);
id_meg_trials = 1:n_trials;

% Atlas
% number of ROIs in the Desikan-Kiliany Atlas
nroi = length(cortex.Atlas(3).Scouts);
% ROI labels
labels = {cortex.Atlas(3).Scouts.Label};
%roi inds 
ind_cortex = [];
ind_roi = {};
for iROI = 1:nroi
  ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices;
  ind_cortex = cat(1, ind_cortex, ind_roi{iROI});
  [~, ind_roi_cortex{iROI}, ~] = intersect(ind_cortex, ind_roi{iROI});
end
nvox = length(ind_cortex);
leadfield = leadfield(:, ind_cortex, :);

%leadfield
% load(sprintf('BF_Patient%s.mat',patientID{id}));
% L = fp_get_lf(inverse);
L = leadfield;
L1 = squeeze(L(:,[inode1 inode2],idir));
L_noise = squeeze(L(:,[rand_nodes],idir));
nmeg = size(L1,1);
id_meg_chan = 1:nmeg;

%project to sensors and add white noise 
for itrial = 1:n_trials
    clear sig whitenoise noise
    sig1 = L1 * xx(:,:,itrial);
    sig1 = sig1 ./ norm(sig1, 'fro');
    noise = L_noise * xx_noise(:,:,itrial);
    noise = noise ./ norm(noise, 'fro');    
    sig = alph * sig1 + (1-alph)* noise;
%     sig=sig1;
    
    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal(:,:,itrial) = sig ./ norm(sig, 'fro');
end

clear x1 x2 xx whitenoise sig X L1 noise

%% megmeg pipeline start
%parameters

filtertype= 'd';
imethod = 'sum';
ndim = 3;
npcs = 2;
fres = 75;

id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(signal,fres,'WELCH',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);

% clear mni_pos label code roi_id u_roi_id csroi
% mni_pos = fp_getMNIpos(patientID{id});
% for ii = 1: ns_org
%     [~,~,roi_id(ii)]=fp_get_mni_anatomy_new(mni_pos(ii,:));
% end
% u_roi_id = sort(unique(roi_id));
% nroi = numel(u_roi_id)-1; %because white voxels are not counted

%construct source filter
if strcmp(filtertype,'e')
    A = squeeze(mkfilt_eloreta_v2(L));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(filtertype,'d')
    A=zeros(nmeg,ndim,nvox,nfreq);
    
    for ifrq = 1:nfreq
        cCS = CS(:,:,ifrq);
        lambda = mean(diag(real(cCS)))/100;
        
        CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
        
        for is=1:nvox %iterate across nodes
            Lloc=squeeze(L(:,is,:));
            A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
        end
    end
    fqA = 1:nfreq; %This filter is frequency specific.
    nfqA = nfreq;
end


clear P
for aroi = 1:nroi
    
    %project to source level
    clear A_ CSv
    A_ = A(:, :,ind_roi_cortex{aroi},:);
    nsroi = size(A_,3);
    A_ = reshape(A_, [nmeg, ndim*nsroi, nfqA]);
    
    for ifq = 1:nfreq
        clear csv
        csv = A_(:,:,fqA(ifq))' * CS(:,:,ifq) * A_(:,:,fqA(ifq));
        n = size(csv,2);
        csv(1:n+1:end) = real(diag(squeeze(csv)));
        CSv(ifq,:,:)=csv;
    end
    
    %zscoring
    clear ZS CSz
    ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 1))))./diag(squeeze(sum(real(CSv), 1)))));
    for ifreq = 1:nfreq
        CSz(ifreq,:, :) = ZS'*squeeze(CSv(ifreq,:, :))*ZS;
    end
    
    %region pca
    clear CSs v v5 in V_ D_
    CSs = squeeze(sum(CSz,1)); %covariance
    [V_, D_] = eig(real(CSs));
    [D_, in] = sort(real(diag(D_)), 'descend');
    % variance explained
%     vx_ = cumsum(D_)./sum(D_);
%     varex{aroi} = vx_;
    
    V{aroi} = V_(:,in(1:npcs)); %nregionvoxels*2 x npcs
    
    
    %     %concatenate filters
    for ifq = 1:nfqA
        P(:, :, aroi,ifq) = A_(:,:,fqA(ifq)) * ZS * real(V{aroi});
    end
end

%apply all filters
CSroi = [];
for ifreq = 1:nfreq
    CSroi(:, :, ifreq) = reshape(P(:,:,:,fqA(ifreq)), nmeg, [])'*CS(:, :, ifreq)*reshape(P(:,:,:,fqA(ifreq)), nmeg, []);
end

%divide by power to obtain coherence
clear Cohroi
for ifreq = 1: nfreq
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end

%integrate across npcs
if strcmp(imethod,'sum')
    
    %sum up coherence across npcs
    ic = 1;
    for iroi = 1:nroi
        jc = 1;
        for jroi = 1:nroi
            true_coh(iroi,jroi,:) = squeeze(sum(sum(Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,:),1),2));
            jc = jc+npcs;
        end
        ic=ic+npcs;
    end
    
    
elseif strcmp(imethod,'mim')
    regu=.000001;
    
    ic = 1;
    for ii=1:nroi
        jc=1;
        for jj= 1: nroi
            for ifq = 1: nfreq
                
                cs_red=[];
                cs_red{1} = Cohroi(ic:ic+npcs-1,ic:ic+npcs-1,ifq);
                cs_red{2} = Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,ifq);
                cs_red{3} = Cohroi(jc:jc+npcs-1,jc:jc+npcs-1,ifq);
                
                caainv=inv(real(cs_red{1})+regu*eye(npcs)*mean(diag(real(cs_red{1}))));
                cab=imag(cs_red{2});
                cbbinv=inv(real(cs_red{3})+regu*eye(npcs)*mean(diag(real(cs_red{3}))));
                X=cab*cbbinv*cab';
                % MIM Ewald Eq. 14
                mim(ii,jj,ifq)=(trace(caainv*X));
                caainvsqrt=sqrtm(caainv);
                Y=caainvsqrt*X*caainvsqrt; %Eq. 23
                [~,s,~]=svd(Y);
                % MIC
                mic(ii,jj,ifq)=sqrt(s(1,1));
%                 mem1{ii,ifq}=(trace(caainv*X));           
            end
            jc = jc+npcs;

        end
        ic=ic+npcs;
    end
else
    error('Unknown imethod')
end


%find rois that belong to inode1 and inode2
clear iroi1 iroi2
for iroi = 1:nroi 
    if any(ind_roi{iroi}==inode1)
        iroi1 = iroi;
    elseif any(ind_roi{iroi}==inode2)
        iroi2 = iroi;
    end
    if exist('iroi1','var')& exist('iroi2','var')
        break
    end
end

load cm17

data_in = squeeze(mean(abs(imag(true_coh(:,iroi1,:),3))));
allplots_cortex_BS(cortex_highres, data_in, [min(data_in) max(data_in)],...
    cm17a,'.', smooth_cortex,['megmeg_sim_brainstorm_meanroi' num2str(iroi1)]);

data_in2 = zeros(size(data_in));
data_in2(iroi1)=1;
allplots_cortex_BS(cortex_highres, data_in2, [min(data_in2) max(data_in2)],...
    cm17a,'.', smooth_cortex,['megmeg_sim_brainstorm_roi' num2str(iroi1)]);

data_in = squeeze(mean(abs(imag(true_coh(:,iroi2,:),3))));
allplots_cortex_BS(cortex_highres, data_in, [min(data_in) max(data_in)],...
    cm17a,'.', smooth_cortex,['megmeg_sim_brainstorm_meanroi' num2str(iroi2)]);

data_in2 = zeros(size(data_in));
data_in2(iroi2)=1;
allplots_cortex_BS(cortex_highres, data_in2, [min(data_in2) max(data_in2)],...
    cm17a,'.', smooth_cortex,['megmeg_sim_brainstorm_roi' num2str(iroi2)]);

% 
% 
% data_in3 = zeros(size(data_in));
% data_in3=1:numel(data_in);
% allplots_cortex_BS(cortex_highres, data_in3, [min(data_in3) max(data_in3)],...
%     cm17a,'.', smooth_cortex,['test' num2str(iroi2)]);
