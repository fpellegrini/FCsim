function pow = fp_project_power(wholeCS, filter)

[nmeg, ns, nfreq] = size(filter); 

pow = nan(ns,nfreq);

for ifreq = 1: nfreq    
    
    cfilter = filter(:,:,ifreq)';
    
    for is = 1:ns    
        pow(is,ifreq) = real(cfilter(is,:) * wholeCS(:,:,ifreq) * cfilter(is,:)');        
    end    
end 