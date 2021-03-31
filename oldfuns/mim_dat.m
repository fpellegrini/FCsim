clear all 

Nepo = 200;
Lepo = 100;
N = Nepo*Lepo;
chan = 10;
lag = 5; 
fres = 40;

s1 = randn(chan, N);
s2 = randn(chan, N);
s2(1, :) = circshift(s1(1, :), lag, 2);

iroi = 1; 
jroi = 2; 

ic = 1+ (iroi-1)*chan;
jc = 1+ (jroi-1)*chan; 
regu=.000001;

%% Test 1 

signal = [s1; s2];
% signal_epo = reshape(signal, 2*chan, Lepo, Nepo);
% CS = fp_tsdata_to_cpsd(signal_epo, fres,'WELCH',1:2*chan, 1:2*chan,1:Nepo, 1:Nepo);
S = tsdata_to_cpsd(signal,fres,'WELCH');
S(:,:,1) = [];
% CS(:,:,1) = [];
CSroi = S;
clear Cohroi

%divide by power to obtain coherence
for ifreq = 1: fres
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end

[coh_true_case1 ind_f_case1] = max(squeeze(abs(imag(Cohroi(1, chan+1, :)))))
ifq = ind_f_case1; 
         
cs_red=[];
cs_red{1} = Cohroi(ic:ic+chan-1,ic:ic+chan-1,ifq); %Caa
cs_red{2} = Cohroi(ic:ic+chan-1,jc:jc+chan-1,ifq); %Cab
cs_red{3} = Cohroi(jc:jc+chan-1,jc:jc+chan-1,ifq); %Cbb

caainv=inv(real(cs_red{1})+regu*eye(chan)*mean(diag(real(cs_red{1}))));
cab=imag(cs_red{2});
cbbinv=inv(real(cs_red{3})+regu*eye(chan)*mean(diag(real(cs_red{3}))));
X=cab*cbbinv*cab';
% MIM Ewald Eq. 14
mim_case1=(trace(caainv*X))
caainvsqrt=sqrtm(caainv);
Y=caainvsqrt*X*caainvsqrt; %Eq. 23
[~,s,~]=svd(Y);
% MIC
mic_case1=sqrt(s(1,1))

%% Test 2
p1 = randn(chan, chan); 
p1 = p1 / norm(p1);
p2 = randn(chan, chan); 
p2 = p2 / norm(p2);

s11 = p1 * s1;
s22 = p2 * s2;
signal = [s11; s22];
% signal_epo = reshape(signal, 2*chan, Lepo, Nepo);
% CS = fp_tsdata_to_cpsd(signal_epo, fres,'WELCH',1:2*chan, 1:2*chan,1:Nepo, 1:Nepo);

S11 = tsdata_to_cpsd(signal,fres,'WELCH');
S11(:,:,1) = [];
% CS(:,:,1) = [];
CSroi = S11;
clear Cohroi

%divide by power to obtain coherence
for ifreq = 1: fres
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end

[coh_true_case2 ind_f_case2] = max(squeeze(abs(imag(Cohroi(1, chan+1, :)))))
ifq = ind_f_case1; 
         
cs_red=[];
cs_red{1} = Cohroi(ic:ic+chan-1,ic:ic+chan-1,ifq); %Caa
cs_red{2} = Cohroi(ic:ic+chan-1,jc:jc+chan-1,ifq); %Cab
cs_red{3} = Cohroi(jc:jc+chan-1,jc:jc+chan-1,ifq); %Cbb

caainv=inv(real(cs_red{1})+regu*eye(chan)*mean(diag(real(cs_red{1}))));
cab=imag(cs_red{2});
cbbinv=inv(real(cs_red{3})+regu*eye(chan)*mean(diag(real(cs_red{3}))));
X=cab*cbbinv*cab';
% MIM Ewald Eq. 14
mim_case2=(trace(caainv*X))
caainvsqrt=sqrtm(caainv);
Y=caainvsqrt*X*caainvsqrt; %Eq. 23
[~,s,~]=svd(Y);
% MIC
mic_case2=sqrt(s(1,1))

%% Test 3

whitenoise = randn(size(signal));
whitenoise = whitenoise ./ norm(whitenoise, 'fro');
signal1 = 0.95 *signal + 0.05 * whitenoise;
signal1 = signal1 ./ norm(signal1, 'fro');
% signal_epo = reshape(signal1, 2*chan, Lepo, Nepo);
% CS = fp_tsdata_to_cpsd(signal_epo, fres,'WELCH',1:2*chan, 1:2*chan,1:Nepo, 1:Nepo);

S2 = tsdata_to_cpsd(signal,fres,'WELCH');
S2(:,:,1) = [];
% CS(:,:,1) = [];
CSroi = S2;
clear Cohroi

%divide by power to obtain coherence
for ifreq = 1: fres
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end

[coh_true_case3 ind_f_case3] = max(squeeze(abs(imag(Cohroi(1, chan+1, :)))))
ifq = ind_f_case1; 
         
cs_red=[];
cs_red{1} = Cohroi(ic:ic+chan-1,ic:ic+chan-1,ifq); %Caa
cs_red{2} = Cohroi(ic:ic+chan-1,jc:jc+chan-1,ifq); %Cab
cs_red{3} = Cohroi(jc:jc+chan-1,jc:jc+chan-1,ifq); %Cbb

caainv=inv(real(cs_red{1})+regu*eye(chan)*mean(diag(real(cs_red{1}))));
cab=imag(cs_red{2});
cbbinv=inv(real(cs_red{3})+regu*eye(chan)*mean(diag(real(cs_red{3}))));
X=cab*cbbinv*cab';
% MIM Ewald Eq. 14
mim_case3=(trace(caainv*X))
caainvsqrt=sqrtm(caainv);
Y=caainvsqrt*X*caainvsqrt; %Eq. 23
[~,s,~]=svd(Y);
% MIC
mic_case3=sqrt(s(1,1))



%% Test 4 

signal4 = [s1; s2];

signal4([ic jc],:) = signal4([ic jc],:)./ norm(signal4([ic jc],:),'fro');
signal4([(ic+1):(jc-1),(jc+1):end],:) = signal4([(ic+1):(jc-1),(jc+1):end],:)./ ...
    norm(signal4([(ic+1):(jc-1),(jc+1):end],:),'fro');
S4 = tsdata_to_cpsd(signal4,fres,'WELCH');
S4(:,:,1) = [];


clear CSroi V
npcs = 5; 
V= zeros(chan*2, npcs*2);
CSs1 = squeeze(sum(S4(1:chan,1:chan,:),3)); %covariance
CSs2 = squeeze(sum(S4(chan+1:end,chan+1:end,:),3)); %covariance
[v1, d1, ~] = eig(real(CSs1));
[v2, d2, ~] = eig(real(CSs2));
[~, in1] = sort(diag(d1), 'descend');
[~, in2] = sort(diag(d2), 'descend');

V(1:chan,1:npcs) = v1(:,in1(1:npcs)); 
V(chan+1:end,npcs+1:end) = v2(:,in2(1:npcs));

CSroi = [];
for ifreq = 1:fres
    CSroi(:, :, ifreq) = V'*S4(:, :, ifreq)*V;
end

clear Cohroi

%divide by power to obtain coherence
for ifreq = 1: fres
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end
   
% ifq = 20;
chan = npcs;
ic = 1+ (iroi-1)*chan;
jc = 1+ (jroi-1)*chan; 

cs_red=[];
cs_red{1} = Cohroi(ic:ic+chan-1,ic:ic+chan-1,ifq); %Caa
cs_red{2} = Cohroi(ic:ic+chan-1,jc:jc+chan-1,ifq); %Cab
cs_red{3} = Cohroi(jc:jc+chan-1,jc:jc+chan-1,ifq); %Cbb

caainv=inv(real(cs_red{1})+regu*eye(chan)*mean(diag(real(cs_red{1}))));
cab=imag(cs_red{2});
cbbinv=inv(real(cs_red{3})+regu*eye(chan)*mean(diag(real(cs_red{3}))));
X=cab*cbbinv*cab';
% MIM Ewald Eq. 14
mim_case4=(trace(caainv*X))
caainvsqrt=sqrtm(caainv);
Y=caainvsqrt*X*caainvsqrt; %Eq. 23
[~,s,~]=svd(Y);
% MIC
mic_case4=sqrt(s(1,1))
