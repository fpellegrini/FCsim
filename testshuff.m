o = COH_o{1};
s = COH_s{1};
dc = abs(o) - abs(s);
fr = 1:200;

figure;
imagesc(log(squeeze(mean(abs(o(fr,:,:)),1))));

figure
imagesc(log(squeeze(mean(abs(s(fr,:,:)),1))));

figure
imagesc(zscore(squeeze(mean(dc(fr,:,:),1))));


po = P_o{1};
ps = P_s{1};
dp = po-ps;

figure
imagesc(po)

figure 
imagesc(ps)

figure 
imagesc(dp)

vo = Pv_o{1};
vs=Pv_s{1};
dv = vo-vs;

figure
imagesc(vo)

figure 
imagesc(vs)

figure 
imagesc(dv)

co = CSv_o{1};
cs = CSv_s{1};
dcs = abs(co) - abs(cs);

figure;
imagesc(log(squeeze(mean(abs(co(fr,:,:)),1))));

figure
imagesc(log(squeeze(mean(abs(cs(fr,:,:)),1))));

figure
imagesc((squeeze(mean(dc(fr,:,:),1))));


