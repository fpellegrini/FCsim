function inds = npcs2inds(npcs)


nroi = length(npcs);
inds = {}; ninds = 0;
for iroi = 1:nroi
    for jroi = (iroi+1):nroi
        inds{ninds+1} = {(iroi-1)*npcs(iroi) + [1:npcs(iroi)], (jroi-1)*npcs(jroi) + [1:npcs(jroi)]};
        inds{ninds+2} = {(jroi-1)*npcs(jroi) + [1:npcs(jroi)], (iroi-1)*npcs(iroi) + [1:npcs(iroi)]};
        ninds = ninds + 2;
    end
end