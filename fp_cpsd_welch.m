function S = fp_cpsd_welch(X_1, X_2,ind_1,ind_2,h,window,noverlap)
%X_1 and X_2 can be the same - then the usual CS is
%calculated. Otherwise X_2 should contain data from another trial than
%X_1. ind_1 and ind_2 contain channels of interest. The output S will be of 
%size length(ind_1) x length(ind_2) x nfreq. 

nfft = 2*(h-1);
n1 = numel(ind_1);
n2 = numel(ind_2);
S = complex(zeros(n1,n2,h));

for ii = 1:n1
    o = ind_1(ii); 
    S(ii,:,:) = transpose(cpsd(X_1(:,o),X_2(:,ind_2),window,noverlap,nfft)); % cross-spectra           
end       

% same as:
% S1 = complex(zeros(n1,n2,h));
% for ii = 1:n1
%     o = ind_1(ii);
%     for jj = 1:n2
%         u = ind_2(jj);
%         S1(ii,jj,:) = cpsd(X_1(:,o),X_2(:,u),window,noverlap,nfft); % cross-spectra 
%     end
% end
