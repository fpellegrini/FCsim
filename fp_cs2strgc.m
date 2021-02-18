function [TRGC, GC, nlags, inds] = fp_cs2strgc(CS, fres, nlags, maxfreq, inds)
% Epoched time series data to spectral conditional time-reversed Granger causality
%
% (C) 2018 Stefan Haufe
%
% This code implements/calls methods presented in the follow papers. Please
% cite the appropriate ones if you use it for a publication.
%
% Haufe, S., Nikulin, V. V., M?ller, K. R., & Nolte, G. (2013). A critical
% assessment of connectivity measures for EEG data: a simulation study.
% Neuroimage, 64, 120-133.
%
% Haufe, S., Nikulin, V. V., & Nolte, G. (2012, March). Alleviating the
% influence of weak data asymmetries on granger-causal analyses.
% In International Conference on Latent Variable Analysis and Signal
% Separation (pp. 25-33). Springer, Berlin, Heidelberg.
%
% Winkler, I., Panknin, D., Bartz, D., M?ller, K. R., & Haufe, S. (2016).
% Validity of time reversal for testing Granger causality. IEEE Transactions
% on Signal Processing, 64(11), 2746-2760.
%
% Barnett, L.; Seth, A.K. Granger causality for state-space models.
% Phys. Rev. E 2015, 91, 040101.
%
% Solo, V. State-space analysis of Granger-Geweke causality measures with
% application to fMRI. Neural Computation 2016, 28, 914?949.
%
% Barnett, L., & Seth, A. K. (2014). The MVGC multivariate Granger causality
% toolbox: a new approach to Granger-causal inference. Journal of
% neuroscience methods, 223, 50-68.


% [nchan, ndat, nepo] = size(data);

if nargin < 2 || isempty(fres)
    fres = ndat;
end
freqs = linspace(0, 1, fres+1);

if nargin < 3 || isempty(nlags)
    nlags = -10;
end

if nargin < 4 || isempty(maxfreq)
    maxfreq = fres+1;
end

freqs = freqs(1:maxfreq);
% z = exp(-i*pi*freqs);
z = exp(-i*pi*(freqs./(fs/2)))';

if nargin < 5 || isempty(inds)
    inds = {}; ninds = 0;
    for ii = 1:nchan
        for ij = (ii+1):nchan
            inds{ninds+1} = {ii, ij};
            inds{ninds+2} = {ij, ii};
            ninds = ninds + 2;
        end
    end
else
    ninds = length(inds);
end


if  nlags < 0
    
    nlagsall = [];
    for ii = 1:min(nchan*(nchan-1)/2, -nlags)
        pe = randperm(nchan);
        [~,~,~, nlags_] = tsdata_to_infocrit(data(pe(1:2), :, :), -nlags, [], 0);
        nlagsall = [nlagsall nlags_];
    end
    nlags = floor(median(nlagsall));
end

clear TRGC GC


% data to autocovariance
%     G = tsdata_to_autocov(data, nlags);
G = cpsd_to_autocov(CS, nlags);

% (time-reversed) GC just between sender and receiver sets

% loop over sender/receiver combinations to compute time-reversed GC
for iind = 1:ninds
    disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] -> [' num2str(inds{iind}{2}) ']'])
    
    subset = [inds{iind}{1} inds{iind}{2}];
    nsubsetvars = length(subset);
    subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
    
    % autocovariance to full forward VAR model
    [A, SIG] = autocov_to_var2(G(subset, subset, :));
    
    % forward VAR model to state space VARMA models
    [eA2, eC2, eK2, eV2, eVy2] = varma2iss(reshape(A, nsubsetvars, []), [], SIG, eye(nsubsetvars));
    
    % backward autocovariance to full backward VAR model
    [AR, SIGR] = autocov_to_var2(permute(G(subset, subset, :), [2 1 3]));
    
    % backward VAR to VARMA
    [eA2R, eC2R, eK2R, eV2R, eVy2R] = varma2iss(reshape(AR, nsubsetvars, []), [], SIGR, eye(nsubsetvars));
    
    % GC and TRGC computation
    GC(:, iind) = iss_SGC(eA2, eC2, eK2, eV2, z, subinds{2}, subinds{1});
    GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, subinds{2}, subinds{1});
    TRGC(:, iind) = GC(:, iind) - GCR';
end