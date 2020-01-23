clear all

%% signal generation
%parameters

load('./processed_bs/bs_results.mat')
smooth_cortex = 0.35;

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2;
% nfreq = 46;
delay = 30; %in samples, is equal to 100 ms
inode1 = 2100; %randi(size(A,2),1);
inode2 = 1000;
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

%leadfield
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);
p = randn(2, 1); 
p = p / norm(p);
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
%     sig = alph * sig1 + (1-alph)* noise;
    sig=sig1;
    
    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal(:,:,itrial) = sig ./ norm(sig, 'fro');
end

signal(end+1,:,:) =  x1;

clear x1 x2 xx whitenoise sig X L1 noise

%% megmeg pipeline start
%parameters

fs = D.fsample;
fres = 75;
frqs = sfreqs(fres, fs);
frqs(frqs>90) = [];
nfreq = numel(frqs);
filtertype= 'e';

id_meg_chan = 1:size(signal,1)-1;
nmeg = numel(id_meg_chan);
id_lfp_chan = size(signal,1);
nlfp = numel(id_lfp_chan);
ndim = size(L,3);

id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(signal,fres,'MT',[id_meg_chan id_lfp_chan], [id_meg_chan id_lfp_chan], id_trials_1, id_trials_2);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);

%construct source filter
if strcmp(filtertype,'e')
    A = squeeze(mkfilt_eloreta_v2(L));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(filtertype,'d')
    ns_org = size(L,2);
    A=zeros(nmeg,ndim,ns_org,nfreq);
    
    for ifrq = 1:nfreq
        cCS = CS(1:end-nlfp,1:end-nlfp,ifrq);
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

A = permute(A,[1 3 2 4]);
cCS = CS(1:(end-nlfp),end-nlfp+1:end,:);

for ifq = 1: nfreq
    for idim =1:ndim
        for ifq = 1:nfreq
            CSv(ifq,idim,:,:) = squeeze(A(:,:,idim,fqA(ifq)))' * cCS(:,:,ifq);
        end

        %get voxel power
        pv(idim,:,:) = fp_project_power_2D(CS(1:nmeg,1:nmeg,:),squeeze(A(:,:,idim,fqA(ifq))))';

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
        clear plfp
        plfp = diag(squeeze(CS(nmeg+1:end,nmeg+1:end,ifreq))); %power of lfp channels
        coh(ifreq, idim,:, :) = squeeze(CSv(ifreq, idim,:, :)) ...
            ./ sqrt(squeeze(pv(idim,ifreq,:))*plfp');
    end
end

% %divide by power to obtain coherence
% clear Cohroi
% for ifreq = 1: nfreq
%     clear pow
%     pow = real(diag(CSroi(:,:,ifreq)));
%     Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
% end

coh1 = squeeze(sum(abs((coh)),2));

a = squeeze(sum(coh1,1)); 

%%

outname = sprintf('true_roi2.nii');
true = zeros(size(a'));
true(inode2)=1;
fp_data2nii(true,sources.pos,[],outname,id)

outname = sprintf('map_megmeg_eloreta.nii');
fp_data2nii(a',sources.pos,[],outname,id)


