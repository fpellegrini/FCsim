

for idim = 1:3 
    for ifq=1:nfreq
    
    u(:,:,idim,ifq) = squeeze(A(:,idim,:))' * CS(:,:,ifq) * squeeze(A(:,idim,:));
    
    end
end

u1 = squeeze(sum(u(:,:,:,filt.iband),4));

plot(abs(imag(u1(:))))

u2 = sum(u1,3);

for aroi = 1:D.nroi 
    for broi = 1:D.nroi
    
        u4(aroi,broi) = sum(sum(abs(imag(u2(D.ind_roi_cortex{aroi},D.ind_roi_cortex{broi})))));
    
    end
end


imagesc(u4)