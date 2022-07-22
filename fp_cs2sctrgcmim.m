function conn = fp_cs2sctrgcmim(CS, fres, nlags, inds, output)

%set parameters 
ninds = length(inds);
freqs = linspace(0, 1, fres+1);
freqs = freqs(1:fres+1);
z = exp(-i*pi*freqs);

if ~isempty(intersect(output, {'MIM', 'MIC', 'COH'}))
    clear COH
    for ifreq = 1:fres+1
        clear pow
        pow = real(diag(CS(:,:,ifreq)));
        COH(:,:,ifreq) = CS(:,:,ifreq)./ sqrt(pow*pow');
    end
end

if ~isempty(intersect(output, {'GC', 'TRGC'}))
    G = cpsd_to_autocov(CS, nlags);
end

for iind = 1:ninds
    if ~isequal(inds{iind}{1}, inds{iind}{2})
        
        %ind configuration
        subset = [inds{iind}{1} inds{iind}{2}];
        nsubsetvars = length(subset);
        subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
        
        if ~isempty(intersect(output, {'MIM', 'MIC'}))
            %MIC and MIM
            [MIC(:, iind) , MIM(:, iind)] =  roi_mim2(COH(subset, subset, :), subinds{1}, subinds{2});
        end
        
        if ~isempty(intersect(output, {'GC', 'TRGC'}))
            % autocovariance to full forward VAR model
            [A, SIG] = autocov_to_var3(G(subset, subset, :));
            
            % forward VAR model to state space VARMA models
            [eA2, eC2, eK2, eV2, eVy2] = varma2iss(reshape(A, nsubsetvars, []), [], SIG, eye(nsubsetvars));
            
            % GC and TRGC computation
            GC(:, iind, 1) = iss_SGC(eA2, eC2, eK2, eV2, z, subinds{2}, subinds{1});
            GC(:, iind, 2) = iss_SGC(eA2, eC2, eK2, eV2, z, subinds{1}, subinds{2});
            
            if ~isempty(intersect(output, {'TRGC'}))
                % backward autocovariance to full backward VAR model
                [AR, SIGR] = autocov_to_var3(permute(G(subset, subset, :), [2 1 3]));
                
                % backward VAR to VARMA
                [eA2R, eC2R, eK2R, eV2R, eVy2R] = varma2iss(reshape(AR, nsubsetvars, []), [], SIGR, eye(nsubsetvars));
                GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, subinds{2}, subinds{1});
                TRGC(:, iind, 1) = GC(:, iind, 1) - GCR';
                GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, subinds{1}, subinds{2});
                TRGC(:, iind, 2) = GC(:, iind, 2) - GCR';
            end
        end
    else
        if ~isempty(intersect(output, {'GC', 'TRGC'}))
            GC(:, iind, 1:2) = 0;
            TRGC(:, iind, 1:2) = 0;
        end
    end
end

if ~isempty(intersect(output, {'COH'}))
    COH = permute(COH, [3 1 2 4]);
end

clear out
for iout = 1:length(output)
    eval(['conn.' output{iout} ' = ' output{iout} ';'])
end
conn.nlags = nlags;
conn.inds = inds;