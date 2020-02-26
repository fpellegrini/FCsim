
clear all

%% signal generation

% 125 sources in the brain
n_trials = 200;
Lepo = 100;
N = n_trials*Lepo;
ns_org = 125;
lag = 5;
fres = 40;
inode_seed = 1;
inode_tar = 63;
filtertype= 'd';
regu=.000001;


s1 = randn(ns_org, N);
s1(inode_tar, :) = circshift(s1(inode_seed, :), lag, 2);


s1([inode_seed inode_tar],:) = s1([inode_seed inode_tar],:)./ norm(s1([inode_seed inode_tar],:),'fro');
s1([(inode_seed+1):(inode_tar-1),(inode_tar+1):end],:) = s1([(inode_seed+1):(inode_tar-1),(inode_tar+1):end],:)./ ...
    norm(s1([(inode_seed+1):(inode_tar-1),(inode_tar+1):end],:),'fro');

signal = reshape(s1, ns_org, Lepo, n_trials);

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2;
%leadfield
load(sprintf('BF_Patient%s.mat',patientID{id}));
L = fp_get_lf(inverse);

ind = randperm(size(L,2),ns_org);
L = L(:,ind,:);
ni = size(L,3);

p = [0.5; 0.5];
p = p/norm(p);

for in = 1:ns_org
    L1 = squeeze(L(:,in,:));
    L_mix(:,in) = L1*p;
end


%project to sensors and add white noise
for itrial = 1:n_trials
    clear sig whitenoise noise
    sig = L_mix' * signal(:,:,itrial);
    sig = sig ./ norm(sig, 'fro');
    
    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal1(:,:,itrial) = sig ./ norm(sig, 'fro');
end

% clear x1 x2 xx whitenoise sig X L1 noise

%% megmeg pipeline start
%parameters

id_meg_chan = 1:size(signal,1);
nmeg = numel(id_meg_chan);

id_trials_1 = 1:n_trials;
id_trials_2 = 1:n_trials;
CS = fp_tsdata_to_cpsd(signal1,fres,'MT',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);

%construct source filter
if strcmp(filtertype,'e')
    A = squeeze(mkfilt_eloreta_v2(L));
    A = permute(A,[1, 3, 2]);
    fqA = ones(1,nfreq);%only one filter for all freqs.
    nfqA = 1;
    
elseif strcmp(filtertype,'d')
    A=zeros(nmeg,ni,ns_org,nfreq);
    
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
elseif strcmp(filtertype,'l')
    
    
end

A2 = reshape(A,nmeg,ni*ns_org,nfreq);

for ifq = 1: nfreq
    CSv(:,:,ifq) = squeeze(A2(:,:,fqA(ifq)))' * CS(:,:,ifq)...
        * squeeze(A2(:,:,fqA(ifq)));
end

%%
CSs = sum(CSv,3); 
[v,d,~] = eig(real(CSs));
[~, in] = sort(diag(d),'descend'); 

V = v(:, in(1:ns_org));

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

CSroi = [];
for ifreq = 1:fres
    CSroi(:, :, ifreq) = V'*CSv(:, :, ifreq)*V;
end

clear Cohroi

%divide by power to obtain coherence
for ifreq = 1: fres
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end
   

%%
% chan = npcs;
% 
% ifq=28;
% cs_red=[];
% cs_red{1} = Cohroi(1:chan,1:chan,ifq); %Caa
% cs_red{2} = Cohroi(1:chan,chan+1:end,ifq); %Cab
% cs_red{3} = Cohroi(chan+1:end,chan+1:end,ifq); %Cbb
chan = 62;
ifq=28;
cs_red=[];
cs_red{1} = Cohroi(1:62,1:62,ifq); %Caa
cs_red{2} = Cohroi(1:62,64:end,ifq); %Cab
cs_red{3} = Cohroi(64:end,64:end,ifq); %Cbb

caainv=inv(real(cs_red{1})+regu*eye(chan)*mean(diag(real(cs_red{1}))));
cab=imag(cs_red{2});
cbbinv=inv(real(cs_red{3})+regu*eye(chan)*mean(diag(real(cs_red{3}))));
X=cab*cbbinv*cab';
% MIM Ewald Eq. 14
mim_case5=(trace(caainv*X))
caainvsqrt=sqrtm(caainv);
Y=caainvsqrt*X*caainvsqrt; %Eq. 23
[~,s,~]=svd(Y);
% MIC
mic_case5=sqrt(s(1,1))
