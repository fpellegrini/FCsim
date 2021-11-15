clear all
DIROUT1=[];
params.iReg=1; %number of interacting voxels in interacting regions 
params.iInt = 1; %number of interacting regions
params.ilag = 2; %lag size
params.isnr = 0.9; %SNR
params.iss = 0.5; %noise mix
params.ip=7; %paramenter configuration 
params.ifilt='e';
params.dimred = 'p';

t=[]; %run time

%%





sig_e = reshape(reshape (signal_roi,ndim,[],200,90),3,68,[]);

%%

for idim = 1:3 
    
   [Pxx_e(idim,:,:),F] = pwelch(squeeze(sig_e(idim,:,:))',fres*2,fres,fres*2,fres);
end 


p_e = squeeze(mean(Pxx_e,1));

plot(F,p_e)

%%

sig_l = reshape(reshape (signal_roi,ndim,[],200,90),3,68,[]);
%

for idim = 1:3 
   
   [Pxx_l(idim,:,:),F] = pwelch(squeeze(sig_l(idim,:,:))',fres*2,fres,fres*2,fres);
end 


p_l = squeeze(mean(Pxx_l,1));

figure
plot(F,p_l)

alpha = find(F>8 & F<12);

%

p_l = squeeze(mean(mean(Pxx_l(:,alpha,:),1),2));
figure; plot(p_l)
%%
p_e = squeeze(mean(mean(Pxx_e(:,alpha,:),1),2));
figure; plot(p_e)