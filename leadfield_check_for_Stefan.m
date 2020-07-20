function leadfield_check_for_Stefan

load('./leadfield_data_for_Stefan.mat')
ivox = 2665;
ntrial = size(X,3);

timeseries = squeeze(X(1,:,:));

%%
ndim = size(L,3);

l = squeeze(L(:,ivox,:));
p = randn(ndim,1);
p = p/norm(p);
lp = l*p;

for itrial = 1:ntrial 
    clear sens whitenoise
    sens =  squeeze(lp* timeseries(:,itrial)');
    whitenoise = randn(size(sens));
    whitenoise = whitenoise ./ norm(whitenoise, 'fro');
    signal(:,:,itrial) = 0.9*sens + 0.1*whitenoise;
end


% %dics
% id_meg_chan = 1:size(signal,1);
% nmeg = numel(id_meg_chan);
% 
% id_trials_1 = 1:ntrial;
% id_trials_2 = 1:ntrial;
% CS = fp_tsdata_to_cpsd(signal,75,'WELCH',[id_meg_chan], [id_meg_chan], id_trials_1, id_trials_2);
% 
% ns_org = size(L,2);
% nfreq = size(CS, 3); 
% ndim = 2; 
% A=zeros(nmeg,ndim,ns_org,nfreq);
% 
% for ifrq = 1:nfreq
%     cCS = CS(:,:,ifrq);
%     lambda = mean(diag(real(cCS)))/100;
%     
%     CSinv=pinv(real(cCS)+lambda * eye(size(cCS)));
%     
%     for is=1:ns_org %iterate across nodes
%         Lloc=squeeze(L(:,is,:));
%         A(:,:,is,ifrq) = (pinv(Lloc'*CSinv*Lloc)*Lloc'*CSinv)'; %create filter
%     end
% end
% A = permute(A,[1 3 2 4]);
% 
% for idim = 1:2 
%     for ifreq = 1:nfreq
%         clear s
%         s = squeeze(A(:,:,idim,ifreq))' * squeeze(CS(:,:,ifreq)) * squeeze(A(:,:,idim,ifreq)); 
%         d_pow(:,ifreq,idim) = diag(s); 
%     end 
% end 
%         
% pow = sum((real(d_pow)),3)'; 


% %eloreta
A = squeeze(mkfilt_eloreta_v2(L));
A = permute(A,[1, 3, 2]);
for idim = 1:2 
    for itrial = 1:size(signal,3)
        clear a s
        a = squeeze(A(:,idim,:)); 
        s = a' * squeeze(signal(:,:,itrial));
        [Pxx(:,:,idim,itrial), f] = pwelch(s', [], [], fs, fs);
    end
end     
pow = sum(sum(Pxx,3),4);

%%
pow(1,:) = []; %cut freq = 0 
alpha = sum(pow(8:12,:),1);
% beta = sum(pow(13:30,:),1);
% gamma = sum(pow(30:end,:),1);

fp_plot_slices_for_Stefan(alpha, pos, 3, './alpha', [0 3] ,mask)