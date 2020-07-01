function npcs = fp_symmetrize_hemispheres(D,npcs,mode1)

partner_rois(1,:) = 1:D.nroi;
partner_rois(2,[1:2:D.nroi-1])=2:2:D.nroi;
partner_rois(2,[2:2:D.nroi]) = 1:2:D.nroi-1;

if strcmp(mode1,'max')
    for iroi = 1:D.nroi
        npcs(iroi) = min(npcs(iroi),npcs(partner_rois(2,iroi)));
    end
    
elseif strcmp(mode1,'percent')
    for iroi = 1:D.nroi
        npcs(iroi) = max(npcs(iroi),npcs(partner_rois(2,iroi)));
    end
end