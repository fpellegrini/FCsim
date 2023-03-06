function [MIM, MIC, iCOH, aCOH] = fp_correct_mim(signal_roi,inds,fres,D,filt,PCA_inds)

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

nit = 50;
[nchan, ~, nepo] = size(signal_roi);

%% true MIM and MIC
output = {'MIM','MIC','COH'};
conn = data2sctrgcmim(signal_roi, fres, 30, 0,0, [], inds, output);
[MIM_t, MIC_t, ~, ~, iCOH_t, aCOH_t] = fp_unwrap_conn(conn,D.nroi,filt,PCA_inds);

%% null distribution

maxfreq = fres+1;
ninds = length(inds);
    
for iit = 1:nit %one iteration takes ~90 sec on my local laptop

    %shuffle trials
    shuf_inds = randperm(nepo);   
    
    clear MIM MIC CS COH
    CS = fp_tsdata_to_cpsd(signal_roi, fres, 'WELCH', 1:nchan, 1:nchan,1:nepo,shuf_inds);
    
    for ifreq = 1:maxfreq
        clear pow
        pow = real(diag(CS(:,:,ifreq)));
        COH(:,:,ifreq) = CS(:,:,ifreq)./ sqrt(pow*pow');
    end
    
    % loop over sender/receiver combinations to compute time-reversed GC
    for iind = 1:ninds
        if ~isequal(inds{iind}{1}, inds{iind}{2})
            %ind configuration
            subset = [inds{iind}{1} inds{iind}{2}];
            subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
            
            %MIC and MIM
            [MIC(:, iind) , MIM(:, iind)] =  roi_mim2(COH(subset, subset, :), subinds{1}, subinds{2});
        end
    end
    
    COH = permute(COH, [3 1 2 4]);   
    
    clear out
    for iout = 1:length(output)
        eval(['conn.' output{iout} ' = ' output{iout} ';'])
    end
    conn.inds = inds;
    
    % extract measures out of the conn struct
    [MIM_s(:,:,iit), MIC_s(:,:,iit), ~, ~, iCOH_s(:,:,iit), aCOH_s(:,:,iit)] = fp_unwrap_conn(conn,D.nroi,filt,PCA_inds);

end

%% Normalization

clear m_ s_
m_ = mean(MIC_s,3);
s_ = std(MIC_s,0,3);
MIC = (MIC_t - m_)./s_;
MIC(isnan(MIC))=0; %fill diagonal with 0 instead of nan

clear m_ s_
m_ = mean(MIM_s,3);
s_ = std(MIM_s,0,3);
MIM = (MIM_t - m_)./s_;
MIM(isnan(MIM))=0;

clear m_ s_
m_ = mean(iCOH_s,3);
s_ = std(iCOH_s,0,3);
iCOH = (iCOH_t - m_)./s_;
iCOH(isnan(iCOH))=0;

clear m_ s_
m_ = mean(aCOH_s,3);
s_ = std(aCOH_s,0,3);
aCOH = (aCOH_t - m_)./s_;
aCOH(isnan(aCOH))=0;

