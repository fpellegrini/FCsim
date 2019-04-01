function pow = fp_project_power(wholeCS, filter)

[nmeg, ns, nfreq] = size(filter); 
nlfp = size(wholeCS,1)-nmeg;

pow = nan(ns,nfreq);

for ifreq = 1: nfreq    
    
    cfilter = cat(2,filter(:,:,ifreq)',eye(ns,nlfp));
    
    for is = 1:ns    
        pow(is,ifreq) = cfilter(is,:) * wholeCS(:,:,ifreq) * cfilter(is,:)';        
    end    
end 