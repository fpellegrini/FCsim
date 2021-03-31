clear all

%% signal generation
%parameters

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};
id = 2;
delay = 30; %in samples, is equal to 100 ms
inode1 = 2100; %randi(size(A,2),1);
inode2 = 1000;
n_rand_nodes = 100;
rand_nodes = randi(2400,1,n_rand_nodes); %generate 100 nodes with random signals
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
    
    %add white noise
    whitenoise = randn(size(sig));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    sig = 0.9*sig + 0.1*whitenoise;
    signal(:,:,itrial) = sig ./ norm(sig, 'fro');
end

clear x1 x2 xx whitenoise sig X L L1 noise

%% megmeg pipeline start
%parameters

filtertype= 'd';
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
    
    %sum up coherence across npcs
    ic = 1;
    for iroi = 1:nroi
        jc = 1;
        for jroi = 1:nroi
            true_coh(iroi,jroi,:) = squeeze(sum(sum(abs(imag(Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,:))),1),2));
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
                cs_red{1} = Cohroi(ic:ic+npcs-1,ic:ic+npcs-1,ifq); %Caa
                cs_red{2} = Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,ifq); %Cab
                cs_red{3} = Cohroi(jc:jc+npcs-1,jc:jc+npcs-1,ifq); %Cbb
                
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

%% plot result
mim=mim.*10^3;
mic=mic.*10^3;
%%

tc1 = mean(true_coh(30,:,:),3);
tc2 = mean(true_coh(18,:,:),3);
x1 = zeros(ns_org,1);
x2 = zeros(ns_org,1);

for ii = 1:ns_org   
    if roi_id(ii)> 0
        x1(ii) = tc1(roi_id(ii));
        x2(ii) = tc2(roi_id(ii));
    end   
end



outname = sprintf('megmeg_sim_%s_to_%d_new.nii',filtertype,inode1);
fp_data2nii(x1,sources.pos,[],outname,id)

outname = sprintf('megmeg_sim_%s_to_%d_new.nii',filtertype,inode2);
fp_data2nii(x2,sources.pos,[],outname,id)


% x_mim = zeros(ns_org,ns_org,nfreq);
% x_mic = zeros(ns_org,ns_org,nfreq);
% 
% tic
% for ifreq = 1:nfreq
%     
%     for ii = 1:ns_org
%         for jj=1:ns_org
%             
%             if roi_id(ii)> 0 && roi_id(jj)>0
%                 x_mim(ii,jj,ifreq) = true_coh(roi_id(ii),roi_id(jj),ifreq);
% %                 x_mim(ii,jj,ifreq) = mim(roi_id(ii),roi_id(jj),ifreq);
% %                 x_mic(ii,jj,ifreq) = mic(roi_id(ii),roi_id(jj),ifreq);
%             end
%             
%         end
%     end
%     
% end
% toc
% 
% a = mean(x_mim(inode1,:,:),3);
% b = mean(x_mim(inode2,:,:),3);

% outname = sprintf('1megmeg_sim_%s_trueandrand_mim_to_%d.nii',filtertype,inode2);
% fp_data2nii(abs(imag(a)),sources.pos,[],outname,id)
% 
% outname = sprintf('1megmeg_sim_%s_trueandrand_mim_to_%d.nii',filtertype,inode1);
% fp_data2nii(abs(imag(b)),sources.pos,[],outname,id)

