function [PERFORMANCE, BASELINE] = fp_get_performance(gt, mic, mim,mean_coh)
%performance dimensions: (mim/mic, pipeline,perfomance measure)
%baseline dimensions: (mim/mix, performance measure)

gt.mic = sum(gt.mic,3);
gt.mim = sum(gt.mim,3);

for ii =1:numel(mim.fixed)
    mic.fixed{ii} = sum(mic.fixed{ii},3);
    mim.fixed{ii} = sum(mim.fixed{ii},3);
    mc.fixed{ii} = sum(mean_coh.fixed{ii},3);
    
    mic.fixed_zs0{ii} = sum(mic.fixed_zs0{ii},3);
    mim.fixed_zs0{ii} = sum(mim.fixed_zs0{ii},3);
    mc.fixed_zs0{ii} = sum(mean_coh.fixed_zs0{ii},3);
    
end
mic.max = sum(mic.max,3);
mim.max = sum(mim.max,3);
mc.max = sum(mean_coh.max,3);
mic.max_zs0 = sum(mic.max_zs0,3);
mim.max_zs0 = sum(mim.max_zs0,3);
mc.max_zs0 = sum(mean_coh_zs0.max,3);
mic.percent = sum(mic.percent,3);
mim.percent = sum(mim.percent,3);
mc.percent = sum(mean_coh.percent,3);
mic.percent_zs0 = sum(mic.percent_zs0,3);
mim.percent_zs0 = sum(mim.percent_zs0,3);
mc.percent_zs0 = sum(mean_coh.percent_zs0,3);
mic.case2 = sum(mic.case2,3);
mim.case2 = sum(mim.case2,3);
mic.baseline = sum(mic.baseline,3);
mim.baseline = sum(mim.baseline,3);


% (1)correlation mim/mic and ground truth
for ii = 1:numel(mic.fixed)
    PERFORMANCE(1,ii,1) = corr(mic.fixed{ii}(:),gt.mic(:));
    PERFORMANCE(2,ii,1) = corr(mim.fixed{ii}(:),gt.mim(:));
    PERFORMANCE(3,ii,1) = corr(mc.fixed{ii}(:),gt.mic(:)); 
    
    PERFORMANCE(1,ii+10,1) = corr(mic.fixed_zs0{ii}(:),gt.mic(:));
    PERFORMANCE(2,ii+10,1) = corr(mim.fixed_zs0{ii}(:),gt.mim(:));
    PERFORMANCE(3,ii+10,1) = corr(mc.fixed_zs0{ii}(:),gt.mic(:)); 
end
PERFORMANCE(1,6,1) = corr(mic.max(:),gt.mic(:));
PERFORMANCE(1,7,1) = corr(mic.percent(:),gt.mic(:));
PERFORMANCE(1,8,1) = corr(mic.case2(:),gt.mic(:));
PERFORMANCE(1,9,1) = corr(mic.max_zs0(:),gt.mic(:));
PERFORMANCE(1,10,1) = corr(mic.percent_zs0(:),gt.mic(:));

PERFORMANCE(2,6,1) = corr(mim.max(:),gt.mim(:));
PERFORMANCE(2,7,1) = corr(mim.percent(:),gt.mim(:));
PERFORMANCE(2,8,1) = corr(mim.case2(:),gt.mim(:));
PERFORMANCE(2,9,1) = corr(mim.max_zs0(:),gt.mim(:));
PERFORMANCE(2,10,1) = corr(mim.percent_zs0(:),gt.mim(:));

PERFORMANCE(3,6,1) = corr(mc.max(:),gt.mic(:));
PERFORMANCE(3,7,1) = corr(mc.percent(:),gt.mic(:));
PERFORMANCE(3,9,1) = corr(mc.max_zs0(:),gt.mic(:));
PERFORMANCE(3,10,1) = corr(mc.percent_zs0(:),gt.mic(:));

%baseline
BASELINE(1,1) = corr(mic.baseline(:),gt.mic(:));
BASELINE(2,1) = corr(mim.baseline(:),gt.mim(:));

