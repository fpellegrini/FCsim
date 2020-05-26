clear all

%% signal generation
%parameters

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2;
% nfreq = 46;
delay = 30; %in samples, is equal to 100 ms
inode1 = 2100; %randi(size(A,2),1);
inode2 = 1000;
isens = randi(125,1);
idir = randi(2,1);
ind = randperm(70);

%load MEG time series and CS
D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
X = D(:,:,:);
id_meg_chan = 1:125;
X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;
 
x1 = squeeze(X(isens,:,:)); %random time meg series
x2 = squeeze(X(isens+1,:,ind)); %second time series with shuffled trials 

xx(1,:,:)= x1;
xx(2,:,:) = x2;
n_trials = size(x1,2);
id_meg_trials = 1:n_trials;

%leadfield
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);
L1 = squeeze(L(:,[inode1 inode2],idir));
nmeg = size(L1,1);
id_meg_chan = 1:nmeg;

for itrial = 1:n_trials
    clear sig whitenoise
    sig = L1 * xx(:,:,itrial);
    sig = sig ./ norm(sig, 'fro');
    %signal on sensor level
    
    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal(:,:,itrial) = sig ./ norm(sig, 'fro');
end

clear x1 x2 xx whitenoise sig X L L1

%% megmeg pipeline start
%parameters

filtertype= 'e';
imethod = 'sum';
ndim = 2;
npcs = 5;
fres = 75;

id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(signal,fres,'WELCH',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);

%leadfield
clear L
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);
ns_org = size(L,2);

clear mni_pos label code roi_id u_roi_id csroi
mni_pos = fp_getMNIpos(patientID{id});
for ii = 1: ns_org
    [~,~,roi_id(ii)]=fp_get_mni_anatomy_new(mni_pos(ii,:));
end
u_roi_id = sort(unique(roi_id));
nroi = numel(u_roi_id)-1; %because white voxels are not counted

%construct source filter
if strcmp(filtertype,'e')
    A = squeeze(mkfilt_eloreta_v2(L));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(filtertype,'d')
    A=zeros(nmeg,ndim,ns_org,nfreq);
    
    for ifrq = 1:nfreq
        cCS = CS(:,:,ifrq);
        lambda = mean(diag(real(cCS)))/100;
        
        CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
        
        for is=1:ns_org %iterate across nodes
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
    A_ = A(:, :,roi_id == aroi,:);
    nsroi = size(A_,3);
    A_ = reshape(A_, [nmeg, ndim*nsroi, nfqA]);
    
    for ifq = 1:nfreq
        CSv(ifq,:,:) = A_(:,:,fqA(ifq))' * CS(:,:,ifq) * A_(:,:,fqA(ifq));
    end
    
    %zscoring
    clear ZS CSz
    ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 1))))./diag(squeeze(sum(real(CSv), 1)))));
    for ifreq = 1:nfreq
        CSz(ifreq,:, :) = ZS'*squeeze(CSv(ifreq,:, :))*ZS;
    end
    
    %region pca
    clear CSs v v5
    CSs = squeeze(sum(CSz,1)); %covariance
    [v, ~, ~] = eig(real(CSs));
    V{aroi} = v(:,1:npcs); %nregionvoxels*2 x npcs
    
    
    %concatenate filters
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
    %             imethod
    
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
    imethod
else
    error('Unknown imethod')
end

%% plot result

x = zeros(ns_org,ns_org,nfreq);

tic
for ifreq = 1:nfreq
    
    for ii = 1:ns_org
        for jj=1:ns_org
            
            if roi_id(ii)> 0 && roi_id(jj)>0
                x(ii,jj,ifreq) = true_coh(roi_id(ii),roi_id(jj),ifreq);
            end
            
        end
    end
    
end
toc

a = mean(abs(imag(x(inode1,:,:))),3);

outname = 'rand_megmeg_sim_e.nii';
fp_data2nii(a,sources.pos,[],outname,id)

