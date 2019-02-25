function S = fp_cpsd_welch(X,ind_1,ind_2,h,window,noverlap)

nfft = 2*(h-1);
n1 = numel(ind_1);
n2 = numel(ind_2);

S = complex(zeros(n1,n2,h));

if n1 == n2 
    
    if any(ind_1==ind_2)
        for i = 1:n1 
            S(i,i,:) = pwelch(X(:,i),window,noverlap,nfft);          % auto-spectra
            for j = i+1:n1 % so we don't compute cross-spectra twice
                S(i,j,:) = cpsd(X(:,i),X(:,j),window,noverlap,nfft); % cross-spectra
            end
        end
    end
    
else 
     for i = 1:n1
        for j = 1:n2
            S(i,j,:) = cpsd(X(:,ind_1(i)),X(:,ind_2(j)),window,noverlap,nfft); % cross-spectra
        end
     end
    
end 
    

  
   
