function [trgc, gc] = fp_data2gc(data,inds)

morder = 30;

[A,SIG,~] = tsdata_to_var(data,morder);

%convert VAR to autocov sequence
G = var_to_autocov(A,SIG);

for iind = 1:ninds
    if ~isequal(inds{iind}{1}, inds{iind}{2})
        
        subset = [inds{iind}{1} inds{iind}{2}];
        nsubsetvars = length(subset);
        subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
        
        % convert autocov sequence into net GC scores
        GC(:,iind) = autocov_to_mvgc(G,subinds{2},subinds{1}) - autocov_to_mvgc(G,subinds{1},subinds{2});
        
        %transpose autocov sequence to obtain autocov seq of time-reversed data
        Grev = permute(G(:,iind),[2 1 3]);
        
        %compute net GC on time-reversed data
        Grev = autocov_to_mvgc(Grev,subinds{2},subinds{1}) - autocov_to_mvgc(Grev,subinds{1},subinds{2});
        
        %TRGC is simply the difference between GC on original and time-reversed
        %data
        TRGC(:,iind) = (GC(:,iind) - GCrev)/2;
    end
end