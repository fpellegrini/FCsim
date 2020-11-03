iroi_seed
iroi_tar

a = mic.max;
apcs = to_save.max.npcs;
u = zeros(68,68); 

for ii = 1:68 
    for jj = 1:68 
        u(ii,jj) = min(apcs(ii),apcs(jj));
    end
end

a_corr = sum(a,3)./u;
subplot(2,2,1); imagesc(a_corr)
title('mic max, sum across freqs, divided by number of pcs')
subplot(2,2,2); imagesc(u)
title('min(npcs roi1, npcs roi2) in max pipeline')
%%
b = mic.percent;
bpcs=to_save.percent.npcs;
up = zeros(68,68); 

for ii = 1:68 
    for jj = 1:68 
        up(ii,jj) = min(bpcs(ii),bpcs(jj));
    end
end

b_corr = sum(b,3)./up;

subplot(2,2,3); imagesc(b_corr)
title('mic percent, sum across freqs, divided by number of pcs')
subplot(2,2,4); imagesc(up)
title('min(npcs roi1, npcs roi2) in percent pipeline')




%%

iteration1.max.npcs=apcs;
iteration1.percent.npcs=bpcs;
iteration1.max.roh = a;
iteration1.percent.roh = b;
%%
iteration2.max.npcs=apcs;
iteration2.percent.npcs=bpcs;
iteration2.max.roh = a;
iteration2.percent.roh = b;
%%
iteration3.max.npcs=apcs;
iteration3.percent.npcs=bpcs;
iteration3.max.roh = a;
iteration3.percent.roh = b;

%%
asum = sum(a,3);
corr(asum(:),u(:))

bsum = sum(b,3); 
corr(bsum(:),up(:))

figure; imagesc(asum)
figure; imagesc(bsum)


