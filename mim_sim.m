clear all
%% signal generation 
delay = 30; 
D = spm_eeg_load(sprintf('redPLFP%s_off', '07'));
X = D(:,:,:);
fres= 75; 
nchan = 5;
npcs = 5;
 
A = zeros(nchan,1024);
B = zeros(nchan,1024-delay+1);

idx = randi(size(X,3),[(nchan*2)-1,1]);

A = (squeeze(X(20,:,idx(1:nchan)))./10^-6)';
B(1,:) = A(1,delay:end); %second time series with time delay
A = A(:,1:end-delay+1);
B(2:end,:) = (squeeze(X(20,delay:end,idx(nchan+1:end)))./10^-6)';
AB = cat(1,A,B);
AB = AB ./ norm(AB, 'fro');

whitenoise = randn(size(AB));
whitenoise = whitenoise ./ norm(whitenoise, 'fro');

sig =  0.95*AB + 0.05 *whitenoise;
signal = sig ./ norm(sig, 'fro');

% whitenoise = randn(size(AB));
% whitenoise = whitenoise ./ norm(whitenoise, 'fro');
% sig = AB;
% sig = sig ./ norm(sig, 'fro');
% sig([2:nchan, nchan+2:nchan*2],:) = whitenoise([2:nchan, nchan+2:nchan*2],:);
% sig([1, nchan+1],:) =  0.9*sig([1, nchan+1],:) + 0.1 *whitenoise([1, nchan+1],:);
% signal = sig ./ norm(sig, 'fro');

CS = fp_tsdata_to_cpsd(signal,fres,'MT',[1:nchan*2], [1:nchan*2], 1, 1);
CS(:,:,[1 47:end])=[];
nfreq = size(CS,3);

% V= zeros(nchan*2, npcs*2);
% CSs1 = squeeze(sum(CS(1:nchan,1:nchan,:),3)); %covariance
% [v1, ~, ~] = eig(real(CSs1));
% V(1:nchan,1:npcs) = v1(:,1:npcs); 
% CSs2 = squeeze(sum(CS(nchan+1:end,nchan+1:end,:),3)); %covariance
% [v2, ~, ~] = eig(real(CSs2));
% V(nchan+1:end,npcs+1:end) = v2(:,1:npcs);
% 
% CSroi = [];
% for ifreq = 1:nfreq
%     CSroi(:, :, ifreq) = V'*CS(:, :, ifreq)*V;
% end

CSroi=CS;


%divide by power to obtain coherence
for ifreq = 1: nfreq
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end



%%
iroi = 1; 
jroi = 2; 
ifq = 10;

ic = 1+ (iroi-1)*npcs;
jc = 1+ (jroi-1)*npcs; 
regu=.000001;
         
cs_red=[];
cs_red{1} = Cohroi(ic:ic+npcs-1,ic:ic+npcs-1,ifq); %Caa
cs_red{2} = Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,ifq); %Cab
cs_red{3} = Cohroi(jc:jc+npcs-1,jc:jc+npcs-1,ifq); %Cbb

caainv=inv(real(cs_red{1})+regu*eye(npcs)*mean(diag(real(cs_red{1}))));
cab=imag(cs_red{2});
cbbinv=inv(real(cs_red{3})+regu*eye(npcs)*mean(diag(real(cs_red{3}))));
X=cab*cbbinv*cab';
% MIM Ewald Eq. 14
mim=(trace(caainv*X));
caainvsqrt=sqrtm(caainv);
Y=caainvsqrt*X*caainvsqrt; %Eq. 23
[~,s,~]=svd(Y);
% MIC
mic=sqrt(s(1,1));


%%
Caa = Cohroi(ic:ic+npcs-1,ic:ic+npcs-1,ifq); %Caa
Cab = Cohroi(ic:ic+npcs-1,jc:jc+npcs-1,ifq); %Cab
Cbb = Cohroi(jc:jc+npcs-1,jc:jc+npcs-1,ifq); %Cbb
Cba = Cohroi(jc:jc+npcs-1,ic:ic+npcs-1,ifq); %Cab
C = [Caa Cab; Cba Cbb];

T = zeros(size(C));
T(1:npcs,1:size(Caa,2))= (real(Caa))^(-0.5);
T(end-size(Cbb,1)+1:end,end-size(Cbb,2)+1:end)= (real(Cbb))^(-0.5);

D = T*C*T;

E = imag(D(1:size(Caa,1), size(Caa,1)+1:end));

[v1,d1] = eig(E*E'); %alpha
[v2,d2] = eig(E'*E); %beta

%norm(diag(d1-d2),'fro') % -> same d

[~, in] = sort(diag(d1), 'descend');
alpha = v1(:,in(1));
beta = v2(:,in(1));

%%

Coh_max = (alpha'* E * beta)/ (norm(alpha)*norm(beta)); %Eq 7; same result as mic(but inverted sign?)
Coh_allcomps = v1'*E* v2; % abs(diag) increases on diag

ediag = diag(E);
e=sum((ediag))
abs(mic)> sum(ediag) %yes
abs(mim)

