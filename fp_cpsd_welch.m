function S = fp_cpsd_welch(X,ind_1,ind_2,h,window,noverlap)

nfft = 2*(h-1);
n1 = numel(ind_1);
n2 = numel(ind_2);

S = complex(zeros(n1,n2,h));

if n1 == n2 
    
    if any(ind_1==ind_2)
        for i = 1:n1
            o = ind_1(i); 
            S(i,i,:) = pwelch(X(:,o),window,noverlap,nfft);          % auto-spectra
            for j = i+1:n1 % so we don't compute cross-spectra twice
                u = ind_1(j);
                S(i,j,:) = cpsd(X(:,o),X(:,u),window,noverlap,nfft); % cross-spectra
            end
        end
    end
    
else 
     for i = 1:n1
         o = ind_1(i); 
        for j = 1:n2
            u = ind_2(j);
            S(i,j,:) = cpsd(X(:,o),X(:,u),window,noverlap,nfft); % cross-spectra
        end
     end
    
end 
    

  
   
