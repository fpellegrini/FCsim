function S = fp_cpsd_welch(X,ind_1,ind_2,h,window,noverlap)

nfft = 2*(h-1);
n1 = numel(ind_1);
n2 = numel(ind_2);
S = complex(zeros(n1,n2,h));

for i = 1:n1
    o = ind_1(i); 
    S(i,:,:) = transpose(cpsd(X(:,o),X(:,ind_2),window,noverlap,nfft)); % cross-spectra           
end       

% same as:
S1 = complex(zeros(n1,n2,h));
for ii = 1:n1
    for jj = 1:4
        S1(ii,jj,:) = cpsd(X(:,ii),X(:,jj),window,noverlap,nfft); % cross-spectra 
    end
end
