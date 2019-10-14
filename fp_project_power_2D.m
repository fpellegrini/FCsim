function pow = fp_project_power_2D(wholeCS, filter)

[nmeg, ns,~] = size(filter); 
nfreq = size(wholeCS,3);
pow = nan(ns,nfreq);

for ifreq = 1: nfreq    
 
    for is = 1:ns    
        pow(is,ifreq) = real(filter(:,is)' * wholeCS(:,:,ifreq) * filter(:,is)); 
    end    
end 