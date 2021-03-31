function [TRGC, GC, nlags, inds] = data2sctrgc_fp(data, fres, nlags, cond, nboot, maxfreq, inds) 
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
% influence of weak data asymmetries on Granger-causal analyses. 
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

rng('shuffle')

[nchan, ndat, nepo] = size(data);

if nargin < 2 || isempty(fres)
  fres = ndat;
end
freqs = linspace(0, 1, fres+1);

if nargin < 3 || isempty(nlags)
  nlags = -10;
end

if nargin < 4 || isempty(cond)
  cond = 0;
end

if nargin < 5 || isempty(nboot)
  nboot = 0;
else
  if nepo < 2
    warning('Not enough epochs to perform bootstrap.')
  end
end

if nargin < 6 || isempty(maxfreq)
  maxfreq = fres+1;
end

freqs = freqs(1:maxfreq);
z = exp(-i*pi*freqs);

if nargin < 7 || isempty(inds)
  inds = {}; ninds = 0;
  for ii = 1:nchan
    for ij = (ii+1):nchan
      inds{ninds+1} = {ii, ij};    
      inds{ninds+2} = {ij, ii};  
      ninds = ninds + 2;
    end
  end 
end

%transfer indeces and data 
tic
%get indices
ia=[];
ib =[];
ca = [];
cb=[];
for ii = 1:ninds
    ia = [ia inds{ii}{1}];%first index 
    ca = [ca numel(inds{ii}{1})]; %count how many indices are at this position 
    ib = [ib inds{ii}{2}];%second index
    cb = [cb numel(inds{ii}{2})];
end

% u_inds: (data) index of all occurring channels (sorted)
% t_ia and t_ib: transferred indeces: where indices occur in u_inds
[u_inds] = union(ia, ib);  

oa = 1; 
ob = 1; 
for ii = 1:ninds 
    t_ia = [];
    t_ib = [];  
    for jj = 1:ca(ii)
        t_ia = [t_ia find(u_inds==ia(oa))];
        oa = oa+1; 
    end
    for jj = 1:cb(ii)   
        t_ib = [t_ib find(u_inds== ib(ob))];
        ob = ob+1; 
    end
    
    t_inds{ii}{1} = t_ia;
    t_inds{ii}{2} = t_ib;
end 
%only keep channels of interest
data = data(u_inds,:,:);
nchan = size(data,1);
clear inds
inds= t_inds;
ninds = length(inds);
toc

if  nlags < 0
  if cond 
    [~,~,~, nlags] = tsdata_to_infocrit(data, -nlags, [], 0);
  else
    nlagsall = [];
    for ii = 1:min(nchan*(nchan-1)/2, -nlags)
      pe = randperm(nchan);
      [~,~,~, nlags_] = tsdata_to_infocrit(data(pe(1:2), :, :), -nlags, [], 0);
      nlagsall = [nlagsall nlags_];
    end
    nlags = floor(median(nlagsall));
  end
end

clear TRGC GC

if abs(nboot) == 0 % no bootstrap

  % data to autocovariance
%     G = tsdata_to_autocov(data, nlags);
  CS = tsdata_to_cpsd(data, fres, 'MT', ndat);
  G = cpsd_to_autocov(CS, nlags);
  clear CS

  if cond 
  % (time-reversed) GC conditioned on all other variables

    % autocovariance to full forward VAR model
    [A, SIG] = autocov_to_var3(G);

    % forward VAR model to state space VARMA models
    [eA2, eC2, eK2, eV2, eVy2] = varma2iss(reshape(A, nchan, []), [], SIG, eye(nchan)); 

    % backward autocovariance to full backward VAR model
    [AR, SIGR] = autocov_to_var3(permute(G, [2 1 3]));

    % backward VAR to VARMA
    [eA2R, eC2R, eK2R, eV2R, eVy2R] = varma2iss(reshape(AR, nchan, []), [], SIGR, eye(nchan));

    %% loop over sender/receiver combinations to compute (time-reversed) GC 
    for iind = 1:ninds   
      if ~isequal(inds{iind}{1}, inds{iind}{2})
        disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] -> [' num2str(inds{iind}{2}) '], conditional'])

        GC(:, iind) = iss_SGC(eA2, eC2, eK2, eV2, z, inds{iind}{2}, inds{iind}{1});
        GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, inds{iind}{2}, inds{iind}{1});
        TRGC(:, iind) = GC(:, iind) - GCR';
      else
        GC(:, iind) = 0;
        TRGC(:, iind) = 0;
      end
    end      
  else
  % (time-reversed) GC just between sender and receiver sets

    % loop over sender/receiver combinations to compute time-reversed GC 
    for iind = 1:ninds  
      if ~isequal(inds{iind}{1}, inds{iind}{2})
        disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] -> [' num2str(inds{iind}{2}) ']'])

        subset = [inds{iind}{1} inds{iind}{2}];
        nsubsetvars = length(subset);
        subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};      

        % autocovariance to full forward VAR model
        [A, SIG] = autocov_to_var3(G(subset, subset, :));

        % forward VAR model to state space VARMA models
        [eA2, eC2, eK2, eV2, eVy2] = varma2iss(reshape(A, nsubsetvars, []), [], SIG, eye(nsubsetvars)); 

        % backward autocovariance to full backward VAR model
        [AR, SIGR] = autocov_to_var3(permute(G(subset, subset, :), [2 1 3]));

        % backward VAR to VARMA
        [eA2R, eC2R, eK2R, eV2R, eVy2R] = varma2iss(reshape(AR, nsubsetvars, []), [], SIGR, eye(nsubsetvars));

        % GC and TRGC computation
        GC(:, iind) = iss_SGC(eA2, eC2, eK2, eV2, z, subinds{2}, subinds{1});
        GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, subinds{2}, subinds{1});
        TRGC(:, iind) = GC(:, iind) - GCR';
%         
%         F = autocov_to_pwcgc(G(subset, subset, :));
%         FR = autocov_to_pwcgc(permute(G(subset, subset, :), [2 1 3]));
      else
        GC(:, iind) = 0;
        TRGC(:, iind) = 0;
      end
    end
  end



else %bootstrap/ permutation 

  % loop over bootstrap samples/ iterations 
  for iboot = 1:abs(nboot)

      if nboot > 0
          bootinds = randi(nepo, nepo, 1);

          % data to autocovariance
          % G = tsdata_to_autocov(data(:, :, bootinds), nlags);
          CS = tsdata_to_cpsd(data(:, :, boot_inds), fres, 'MT', ndat); 
          
      else           
          % shuffle trials 
          id_trials1 = 1:nepo;
          id_trials2 = randperm(nepo);   
                      
          CS = fp_tsdata_to_cpsd(data, fres, 'MT',1:nchan, 1:nchan, id_trials1, id_trials2, ndat);
      end

    G = cpsd_to_autocov(CS, nlags);

    if cond 
    % (time-reversed) GC conditioned on all other variables

      % autocovariance to full forward VAR model
      [A, SIG] = autocov_to_var3(G);

      % forward VAR model to state space VARMA models
      [eA, eC, eK, eV, eVy] = varma2iss(reshape(A, nchan, []), [], SIG, eye(nchan)); 

      % backward autocovariance to full backward VAR model
      [AR, SIGR] = autocov_to_var3(permute(G, [2 1 3]));

      % backward VAR to VARMA
      [eAR, eCR, eKR, eVR, eVyR] = varma2iss(reshape(AR, nchan, []), [], SIGR, eye(nchan));

      %% loop over sender/receiver combinations to compute (time-reversed) GC 
      for iind = 1:ninds
        if ~isequal(inds{iind}{1}, inds{iind}{2})
%           i3 = setdiff(setdiff(1:nchan, inds{iind}{1}), inds{iind}{2});
          disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] -> [' num2str(inds{iind}{2}) '], conditional'])

          GC(:, iind, iboot) = iss_SGC(eA, eC, eK, eV, z, inds{iind}{2}, inds{iind}{1});
          GCR = iss_SGC(eAR, eCR, eKR, eVR, z, inds{iind}{2}, inds{iind}{1});
          TRGC(:, iind, iboot) = GC(:, iind, iboot) - GCR';
        else
          GC(:, iind, iboot) = 0;
          TRGC(:, iind, iboot) = 0;
        end
      end
    else
    % (time-reversed) GC just between sender and receiver sets

      % loop over sender/receiver combinations to compute time-reversed GC 
      for iind = 1:ninds 
        if ~isequal(inds{iind}{1}, inds{iind}{2})
          disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] -> [' num2str(inds{iind}{2}) ']'])

          subset = [inds{iind}{1} inds{iind}{2}];
          nsubsetvars = length(subset);
          subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};      

          % autocovariance to full forward VAR model
          [A, SIG] = autocov_to_var3(G(subset, subset, :));

          % forward VAR model to state space VARMA models
          [eA, eC, eK, eV, eVy] = varma2iss(reshape(A, nsubsetvars, []), [], SIG, eye(nsubsetvars)); 

          % backward autocovariance to full backward VAR model
          [AR, SIGR] = autocov_to_var3(permute(G(subset, subset, :), [2 1 3]));

          % backward VAR to VARMA
          [eAR, eCR, eKR, eVR, eVyR] = varma2iss(reshape(AR, nsubsetvars, []), [], SIGR, eye(nsubsetvars));

          % GC and TRGC computation
          GC(:, iind, iboot) = iss_SGC(eA, eC, eK, eV, z, subinds{2}, subinds{1});
          GCR = iss_SGC(eAR, eCR, eKR, eVR, z, subinds{2}, subinds{1});
          TRGC(:, iind, iboot) = GC(:, iind, iboot) - GCR';
        else
          GC(:, iind, iboot) = 0;
          TRGC(:, iind, iboot) = 0;
        end
      end        
    end       
  end    
end
  