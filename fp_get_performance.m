function [PERFORMANCE, BASELINE] = fp_get_performance(gt, mic, mim, params)
%performance dimensions: (mim/mic, pipeline,perfomance measure)
%baseline dimensions: (mim/mix, performance measure)

% gt_save = gt;
% mic_save = mic;
% mim_save = mim;

gt.mic = sum(gt.mic,3);
gt.mim = sum(gt.mim,3);

for ii =1:numel(mim.fixed)
    mic.fixed{ii} = sum(mic.fixed{ii},3);
    mim.fixed{ii} = sum(mim.fixed{ii},3);
end
mic.max = sum(mic.max,3);
mim.max = sum(mim.max,3);
mic.percent = sum(mic.percent,3);
mim.percent = sum(mim.percent,3);
mic.case2 = sum(mic.case2,3);
mim.case2 = sum(mim.case2,3);
mic.baseline = sum(mic.baseline,3);
mim.baseline = sum(mim.baseline,3);


% (1)correlation mim/mic and ground truth
for ii = 1:numel(mic.fixed)
    PERFORMANCE(1,ii,1) = corr(mic.fixed{ii}(:),gt.mic(:));
    PERFORMANCE(2,ii,1) = corr(mim.fixed{ii}(:),gt.mim(:));
end
PERFORMANCE(1,6,1) = corr(mic.max(:),gt.mic(:));
PERFORMANCE(1,7,1) = corr(mic.percent(:),gt.mic(:));
PERFORMANCE(1,8,1) = corr(mic.case2(:),gt.mic(:));
PERFORMANCE(2,6,1) = corr(mim.max(:),gt.mim(:));
PERFORMANCE(2,7,1) = corr(mim.percent(:),gt.mim(:));
PERFORMANCE(2,8,1) = corr(mim.case2(:),gt.mim(:));

%baseline
BASELINE(1,1) = corr(mic.baseline(:),gt.mic(:));
BASELINE(2,1) = corr(mim.baseline(:),gt.mim(:));


% (2) correlation maxima of mim/mic and ground truth
for ii = 1:numel(mim.fixed)
    %mic
    clear m_max
    m_max = fp_get_nmaxima(mic.fixed{ii},params.iInt*2);
    PERFORMANCE(1,ii,2) = corr(m_max(:),gt.mic(:));
    
    %mim
    clear m_max
    m_max = fp_get_nmaxima(mim.fixed{ii},params.iInt*2);
    PERFORMANCE(2,ii,2) = corr(m_max(:),gt.mim(:));
end

%max pipeline
%mic
clear m_max
m_max = fp_get_nmaxima(mic.max,params.iInt*2);
PERFORMANCE(1,6,2) = corr(m_max(:),gt.mic(:));
%mim
clear m_max
m_max = fp_get_nmaxima(mim.max,params.iInt*2);
PERFORMANCE(2,6,2) = corr(m_max(:),gt.mim(:));

%percent pipeline
%mic
clear m_max
m_max = fp_get_nmaxima(mic.percent,params.iInt*2);
PERFORMANCE(1,7,2) = corr(m_max(:),gt.mic(:));
%mim
clear m_max
m_max = fp_get_nmaxima(mim.percent,params.iInt*2);
PERFORMANCE(2,7,2) = corr(m_max(:),gt.mim(:));

%case2 pipeline
%mic
clear m_max
m_max = fp_get_nmaxima(mic.case2,params.iInt*2);
PERFORMANCE(1,8,2) = corr(m_max(:),gt.mic(:));
%mim
clear m_max
m_max = fp_get_nmaxima(mim.case2,params.iInt*2);
PERFORMANCE(2,8,2) = corr(m_max(:),gt.mim(:));

%baseline
%mic
clear m_max
m_max = fp_get_nmaxima(mic.baseline,params.iInt*2);
BASELINE(1,2) = corr(m_max(:),gt.mic(:));
%mim
clear m_max
m_max = fp_get_nmaxima(mim.baseline,params.iInt*2);
BASELINE(2,2) = corr(m_max(:),gt.mim(:));