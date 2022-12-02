function fp_simulation_stats2

DIRIN = '/home/bbci/data/haufe/Franziska/data/mim_sim5_stats/';
load([DIRIN 'signal.mat'])

iit = str2num(getenv('SGE_TASK_ID'));
rng(iit)

%% calculate true MIM and do statistics 
tic
nit = 1000;
[nchan, ~, nepo] = size(signal_roi);
maxfreq = fres+1;

% calculate indices
clear inds PCA_inds
[inds, PCA_inds] = fp_npcs2inds(npcs);
ninds = length(inds);
   
for iit = 1:nit %one iteration takes ~90 sec on my local laptop

    %shuffle trials
    shuf_inds = randperm(nepo);   
    
    clear MIM2 CS COH2
    CS = fp_tsdata_to_cpsd(signal_roi, fres, 'WELCH', 1:nchan, 1:nchan,1:nepo,shuf_inds);
    
    for ifreq = 1:maxfreq
        clear pow
        pow = real(diag(CS(:,:,ifreq)));
        COH2(:,:,ifreq) = CS(:,:,ifreq)./ sqrt(pow*pow');
    end
    
    % loop over sender/receiver combinations to compute time-reversed GC
    for iind = 1:ninds
        if ~isequal(inds{iind}{1}, inds{iind}{2})
            %ind configuration
            subset = [inds{iind}{1} inds{iind}{2}];
            subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
            
            %MIC and MIM
            [~ , MIM2(:, iind)] =  roi_mim2(COH2(subset, subset, :), subinds{1}, subinds{2});
        end
    end        
      
    % extract measures out of the conn struct
    clear conn
    conn.MIM = MIM2;
    conn.inds = inds;  
    [MIM_s(:,:,iit), ~, ~, ~, ~, ~] = fp_unwrap_conn(conn,D.nroi,filt,PCA_inds);

end

% for iroi = 1:D.nroi 
%     for jroi = 1:D.nroi
%         MIM_p(iroi,jroi) = sum(squeeze(MIM_s(iroi,jroi,:))>MIM_t(iroi,jroi))/nit;
%     end
% end
t = toc; 
%% save 
outname1 = [DIROUT 'result_' num2str(iit) '.mat'];
save(outname1,'MIM_s','t','-v7.3')




