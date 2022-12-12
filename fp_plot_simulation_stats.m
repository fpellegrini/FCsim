DIRIN = '~/Dropbox/Franziska/Data_MEG_Project/simulation_stats/2/';

load([DIRIN 'signal.mat'])

mim_s = [];
for ii = 1:100 
    clear MIM_s
    load([DIRIN 'result_' num2str(ii) '.mat'])
    mim_s = cat(3,mim_s,MIM_s);
end


for iroi = 1:D.nroi 
    for jroi = 1:D.nroi
        MIM_p(iroi,jroi) = sum(squeeze(mim_s(iroi,jroi,:))>MIM_t(iroi,jroi))/size(mim_s,3);
    end
end

imagesc(-log10(MIM_p))

%%

MIM_p(MIM_p==0) = 0.0001;
for ii = [11 49]
    data = squeeze(-log10(MIM_p(:,ii)));
    data(ii)=-1;
    allplots_cortex_BS(cortex_highres,data, [-max(data) max(data)], cm17 ,'-log10(p)', 0.3,[])
end