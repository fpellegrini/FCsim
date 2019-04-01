function S = cpsd_welch(X,n,h,window,noverlap)

keyboard
nfft = 2*(h-1);

S = complex(zeros(n,n,h));

for i = 1:n %iterate over channels 
    S(i,i,:) = pwelch(X(:,i),window,noverlap,nfft);          % auto-spectra
    for j = i+1:n % so we don't compute cross-spectra twice
        S(i,j,:) = cpsd(X(:,i),X(:,j),window,noverlap,nfft); % cross-spectra
    end
end