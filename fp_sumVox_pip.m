function [MIM_, MIC_] = fp_sumVox_pip(signal_roi, nroi, ndim, nvoxroi, fres, filt)

output = {'MIM','MIC'};

% For memory purposes, we need to calculate this in a
% region-by-region way
for oroi = 1:nroi-1
    for uroi = oroi+1:nroi
        clear data npcs mic mim conn inds
        
        npcs = repmat(ndim,1,nvoxroi(oroi)+nvoxroi(uroi));
        data = cat(1,signal_roi{oroi}, signal_roi{uroi});
        
        [inds, ~] = fp_npcs2inds(npcs);
        
        conn = data2sctrgcmim(data, fres, 20, 0,0, [], inds, output);
        
        % extract measures out of the conn struct
        [mim, mic, ~,~,~] = fp_unwrap_conn(conn,nvoxroi(oroi)+nvoxroi(uroi),filt,[]);
        
        MIM_(oroi,uroi) = squeeze(mean(mean(mim(1:nvoxroi(oroi),nvoxroi(oroi):end),1),2));
        MIM_(uroi,oroi) = MIM_(oroi,uroi);
        MIC_(oroi,uroi) = squeeze(mean(mean(mic(1:nvoxroi(oroi),nvoxroi(oroi):end),1),2));
        MIC_(uroi,oroi) = MIC_(oroi,uroi);
        
    end
end