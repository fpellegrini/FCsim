function S = fp_cpsd_mt(X1,X2,ind_1,ind_2,h,window,noverlap,nchunks,taparray)

nfft = 2*(h-1);

n1 = numel(ind_1);
n2 = numel(ind_2);
S = complex(zeros(h,n1,n2));

winstep = window-noverlap;
ntapers = size(taparray,2);

% compute tapered periodogram with FFT

for k = 1:nchunks
    
    XSEG1 = X1((1:window) + (k-1)*winstep,:);
    XSEG2 = X2((1:window) + (k-1)*winstep,:);
    
    % compute periodogram
    P1 = fft(taparray.*permute(XSEG1(:,:,ones(1,ntapers)),[1 3 2]),nfft);
    P1 = P1(1:h,:,:);
    P2 = fft(taparray.*permute(XSEG2(:,:,ones(1,ntapers)),[1 3 2]),nfft);
    P2 = P2(1:h,:,:);
    
    % now make cross-products of them to fill cross-spectrum matrix
    
    for ii = 1:n1
        o = ind_1(ii);        
        S(:,ii,:) = S(:,ii,:) + mean(P1(:,:,o) .* conj(P2(:,:,ind_2)),2);
    end
    
end

S = S/nchunks;






