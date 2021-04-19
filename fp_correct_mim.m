function [MIM, MIC, iCOH, aCOH] = fp_correct_mim(signal_roi,inds,fres,D,filt,PCA_inds)

nit = 50;
[~,~, n_trials] = size(signal_roi);

%% true MIM and MIC
output = {'MIM','MIC','COH'};
conn = data2sctrgcmim(signal_roi, fres, 20, 0,0, [], inds, output);
% extract measures out of the conn struct
[MIM_t, MIC_t, ~, iCOH_t, aCOH_t] = fp_unwrap_conn(conn,D.nroi,filt,PCA_inds);

%% null distribution

for iit = 1:nit
    
    %shuffle trials
    shuf_inds = randperm(n_trials);
    sig_s = signal_roi(:,:,shuf_inds);
    
    %calculate MIM and MIC
    clear conn
    output = {'MIM','MIC','COH'};
    conn = data2sctrgcmim(sig_s, fres, 20, 0,0, [], inds, output,0);        
    % extract measures out of the conn struct
    [MIM_s(:,:,iit), MIC_s(:,:,iit), ~, iCOH_s(:,:,iit), aCOH_s(:,:,iit)] = fp_unwrap_conn(conn,D.nroi,filt,PCA_inds);    
end

%% Normalization

clear m_ s_
m_ = mean(MIC_s,3);
s_ = std(MIC_s,0,3);
MIC = (MIC_t - m_)./s_;

clear m_ s_
m_ = mean(MIM_s,3);
s_ = std(MIM_s,0,3);
MIM = (MIM_t - m_)./s_;

clear m_ s_
m_ = mean(iCOH_s,3);
s_ = std(iCOH_s,0,3);
iCOH = (iCOH_t - m_)./s_;

clear m_ s_
m_ = mean(aCOH_s,3);
s_ = std(aCOH_s,0,3);
aCOH = (aCOH_t - m_)./s_;

