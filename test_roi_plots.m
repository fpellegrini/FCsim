
x = zeros(ns_org,ns_org,nfreq);

tic
for ifreq = 1:nfreq
    
    for ii = 1:ns_org
        for jj=1:ns_org
            
            if roi_id(ii)> 0 && roi_id(jj)>0
                x(ii,jj,ifreq) = true_coh(roi_id(ii),roi_id(jj),ifreq);
            end
            
        end
    end
    
end
toc


a = mean(abs(imag(x(1000,:,:))),3);
%%

b = zeros(ns_org,1);
for kk = 1:ns_org
    if roi_id(kk)>0
        b(kk) = roi_id(kk); 
    end
end

c = zeros(ns_org,1);
for ll = 1:ns_org
    if roi_id(ll)==0
    c(ll) = 0.8; 
    end
end

inode2
d = zeros(ns_org,1);
for mm = 1:ns_org
    if roi_id(mm)==roi_id(inode2)
        
        d(mm)=0.5;
    end
end

d(inode2)=0.8;


outname = '2testroi.nii';
fp_data2nii(d,sources.pos,[],outname,id)