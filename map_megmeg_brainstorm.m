clear all

filtertype= 'l';
%% signal generation
%parameters

rng(1)
load('./processed_bs/bs_results.mat')
smooth_cortex = 0.35;

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2;
% nfreq = 46;
delay = 30; %in samples, is equal to 100 ms
inode_seed = 10; %randi(size(A,2),1);
inode_tar = 1000;
inodes = [inode_seed inode_tar];
% n_rand_nodes = 100;
% rand_nodes = randi(950,1,n_rand_nodes); %generate 100 nodes with random signals
isens = 5; %randi(125,1); %5 or 53
% alph = 0.5;

%load MEG time series and CS
D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
X = D(:,:,:);
id_meg_chan = 1:125;
X(id_meg_chan,:,:)= X(id_meg_chan,:,:)./10^-6;
n_trials = size(X,3);
 
%correlated nodes 
x1 = squeeze(X(isens,:,:)); %random time meg series
x2 = x1(delay:end,:); %second time series with time delay
x1 = x1(1:end-delay+1,:);
xx(1,:,:)= x1;
xx(2,:,:) = x2;

%uncorrelated nodes
% for ii =1:n_rand_nodes
%     ind = randperm(70);
%     xx_noise(ii,:,:) = squeeze(X(isens,delay:end,ind));
% end

n_trials = size(x1,2);
id_meg_trials = 1:n_trials;

%leadfield
L = leadfield;
ni = size(L,3);

for in = 1:2 % 2 nodes
    clear p L1
    rng(1)
    p = randn(ni, 1); 
    p = p / norm(p);
    L1 = squeeze(L(:,inodes(in),:));
    % L_noise = squeeze(L(:,[rand_nodes],:));
    L_mix(:,in) = L1*p;
end

%project to sensors and add white noise 
for itrial = 1:n_trials
    clear sig sig1 whitenoise noise
    sig1 = L_mix * xx(:,:,itrial);
    sig1 = sig1 ./ norm(sig1, 'fro');
%     noise = L_noise * xx_noise(:,:,itrial);
%     noise = noise ./ norm(noise, 'fro');    
%     sig = alph * sig1 + (1-alph)* noise;
    sig=sig1;

    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal(:,:,itrial) = sig ./ norm(sig, 'fro');
end

clear x1 x2 xx whitenoise sig X L1 noise

%% megmeg pipeline start
%parameters

fs = D.fsample;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90) = [];
nfreq = numel(frqs);

id_meg_chan = 1:size(signal,1);
nmeg = numel(id_meg_chan);

id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(signal,fres,'MT',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);
ns = size(L,2);

for is = 1:ns 
    clear u s L2
    L2 = squeeze(L(:,is,:));
    [u, s, v] = svd(squeeze(L2));
    L3(:,is,:) = u(:,1:2)*s(1:2,1:2);
end
clear L 
L = L3;
clear L2 L3
ndim = size(L,3);

%construct source filter
if strcmp(filtertype,'e')
    A = squeeze(mkfilt_eloreta_v2(L));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    A = permute(A,[1 3 2 4]);
    
elseif strcmp(filtertype,'d')
    ns_org = size(L,2);
    A=zeros(nmeg,ndim,ns_org,nfreq);
    
    for ifrq = 1:nfreq
        clear cCS
        cCS = CS(:,:,ifrq);
%         lambda = mean(diag(real(cCS)))/100;  
%         CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
        reg = 0.05*trace(cCS)/length(cCS);
        Cr = cCS + reg*eye(size(cCS,1));
        CSinv=pinv(real(Cr));
               
        for is=1:ns_org %iterate across nodes
            Lloc=squeeze(L(:,is,:));
            A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
        end
    end
    fqA = 1:nfreq; %This filter is frequency specific.
    nfqA = nfreq;
    A = permute(A,[1 3 2 4]);
    
elseif strcmp(filtertype,'l')
    
    cCS = sum(CS,3);
    reg = 0.05*trace(cCS)/length(cCS);
    Cr = cCS + reg*eye(size(cCS,1));

    [~, A] = lcmv_meg(Cr, L, struct('alpha', 0, 'onedim', 0));
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
end


for ifq = 1: nfreq
    for idim =1:ndim
        for ifq = 1:nfreq
            CSv(ifq,idim,:) = squeeze(A(:,:,idim,fqA(ifq)))' * CS(:,:,ifq)...
                * squeeze(A(:,inode_seed,idim,fqA(ifq)));
        end

        %get voxel power
        pv(idim,:,:) = fp_project_power_2D(CS,squeeze(A(:,:,idim,fqA(ifq))))';

    end
end
    
%     %zscoring
%     clear ZS CSz
%     ZS = diag(sqrt(mean(diag(squeeze(sum(real(CSv), 1))))./diag(squeeze(sum(real(CSv), 1)))));
%     for ifreq = 1:nfreq
%         CSz(ifreq,:, :) = ZS'*squeeze(CSv(ifreq,:, :))*ZS;
%     end

coh = CSv;
for idim = 1:ndim
    for ifreq = 1:nfreq
        coh(ifreq, idim,:, :) = squeeze(CSv(ifreq, idim,:)) ...
            ./ sqrt(squeeze(pv(idim,ifreq,:))*squeeze(pv(idim,ifreq,inode_seed))');
    end
end

coh1 = squeeze(sum(abs(imag(coh)),2)); % sum across dims

a = squeeze(sum(coh1,1)); % sum across freqs

%%

load cm17
pos = cortex.Vertices;

d = eucl(pos,pos(inode_tar,:));
xx = exp(-10^-1.5*d);

data_in=xx;
allplots_cortex_BS(cortex, data_in, [min(data_in) max(data_in)],...
    cm17a,'.', smooth_cortex,['true_target_brainstorm']);
clear data_in

data_in = a';
allplots_cortex_BS(cortex, data_in, [min(data_in) max(data_in)],...
    cm17a,'.', smooth_cortex,['map_megmeg_brainstorm_lcmv_29jan1']);


