clear all 

Nepo = 50;
Lepo = 100;
N = Nepo*Lepo;
chan = 64;
lag = 5;
fres = 40;

% rng(500)
s1 = randn(chan, N);
% rng(2)
s2 = randn(chan, N);
s2(1, :) = circshift(s1(1, :), lag, 2);

%% Test 1 

signal = [s1; s2];

signal_epo = reshape(signal, 2*chan, Lepo, Nepo);


CS = fp_tsdata_to_cpsd(signal_epo, fres,'MT',1:2*chan, 1:2*chan,1:Nepo, 1:Nepo);
S = tsdata_to_cpsd(signal,fres,'MT');

S(:,:,1) = [];
CS(:,:,1) = [];
CSroi = S;

%% Test 2
p1 = randn(chan, chan); 
p1 = p1 / norm(p1);

p2 = randn(chan, chan); 
p2 = p2 / norm(p2);

s11 = p1 * s1;
s22 = p2 * s2;

signal = [s11; s22];

signal_epo = reshape(signal, 2*chan, Lepo, Nepo);


CS = fp_tsdata_to_cpsd(signal_epo, fres,'MT',1:2*chan, 1:2*chan,1:Nepo, 1:Nepo);
S = tsdata_to_cpsd(signal,fres,'MT');

S(:,:,1) = [];
CS(:,:,1) = [];
CSroi = S;

%% Test 3

p1 = randn(chan, chan); 
p1 = p1 / norm(p1);

p2 = randn(chan, chan); 
p2 = p2 / norm(p2);

s11 = p1 * s1;
s22 = p2 * s2;

signal = [s11; s22];

whitenoise = randn(size(signal));
whitenoise = whitenoise ./ norm(whitenoise, 'fro');

signal = 0.95 *signal + 0.05 * whitenoise;
signal = signal ./ norm(signal, 'fro');

signal_epo = reshape(signal, 2*chan, Lepo, Nepo);


CS = fp_tsdata_to_cpsd(signal_epo, fres,'MT',1:2*chan, 1:2*chan,1:Nepo, 1:Nepo);
S = tsdata_to_cpsd(signal,fres,'MT');

S(:,:,1) = [];
CS(:,:,1) = [];
CSroi = S;

%% Test 4 

clear CSroi V
npcs = 5; 
V= zeros(chan*2, npcs*2);
CSs1 = squeeze(sum(S(1:chan,1:chan,:),3)); %covariance
CSs2 = squeeze(sum(S(chan+1:end,chan+1:end,:),3)); %covariance
[v1, d1, ~] = eig(real(CSs1));
[v2, d2, ~] = eig(real(CSs2));
[~, in1] = sort(diag(d1), 'descend');
[~, in2] = sort(diag(d2), 'descend');

V(1:chan,1:npcs) = v1(:,in1(1:npcs)); 
V(chan+1:end,npcs+1:end) = v2(:,in2(1:npcs));


CSroi = [];
for ifreq = 1:fres
    CSroi(:, :, ifreq) = V'*CS(:, :, ifreq)*V;
end

chan = npcs;

%%
% s = abs(imag(S(:,:,5)));
% cs = abs(imag(CS(:,:,5)));
% 
% imagesc(s)
% figure
% imagesc(cs)
% 
% all = reshape(abs(imag(CS)),[],40);
% imagesc(all)

%%
clear Cohroi

%divide by power to obtain coherence
for ifreq = 1: fres
    clear pow
    pow = real(diag(CSroi(:,:,ifreq)));
    Cohroi(:,:,ifreq) = CSroi(:,:,ifreq)./ sqrt(pow*pow');
end

a=abs(imag(Cohroi(:,:,5)));

%%
iroi = 1; 
jroi = 2; 
ifq = 5;

ic = 1+ (iroi-1)*chan;
jc = 1+ (jroi-1)*chan; 
regu=.000001;
         
cs_red=[];
cs_red{1} = Cohroi(ic:ic+chan-1,ic:ic+chan-1,ifq); %Caa
cs_red{2} = Cohroi(ic:ic+chan-1,jc:jc+chan-1,ifq); %Cab
cs_red{3} = Cohroi(jc:jc+chan-1,jc:jc+chan-1,ifq); %Cbb

caainv=inv(real(cs_red{1})+regu*eye(chan)*mean(diag(real(cs_red{1}))));
cab=imag(cs_red{2});
cbbinv=inv(real(cs_red{3})+regu*eye(chan)*mean(diag(real(cs_red{3}))));
X=cab*cbbinv*cab';
% MIM Ewald Eq. 14
mim=(trace(caainv*X));
caainvsqrt=sqrtm(caainv);
Y=caainvsqrt*X*caainvsqrt; %Eq. 23
[~,s,~]=svd(Y);
% MIC
mic=sqrt(s(1,1))