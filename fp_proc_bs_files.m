
saveOutputData = 1;   % 1: To save preprocessd data and preprocessing details  0: Not save

settings.smooth_cortex = 1;

%% Loading EEGLAB

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

study_name = 'Attempt1';

bs_folder = '~/brainstorm_db/';

result_folder = './processed_bs/';
if ~exist(result_folder); mkdir(result_folder); end

% load the cortical surface for which the leadfield has been computed. I
% used a downsampled version with only 2000 nodes, and I used a surface
% My idea is to use an ROI based analyses, so the intial resolution should
% not be crucial. But one could go up to about 5000 nodes.
cortex = load([bs_folder '/Attempt1/anat/@default_subject/tess_cortex_mid_low_1000V.mat']);
nvox = size(cortex.Vertices, 1);

% load a very high resolution version of the same surface for plotting
cortex_highres = load([bs_folder '/Attempt1/anat/@default_subject/tess_cortex_mid_low_5000V.mat']);

% lowres version
cortex_lowres = load([bs_folder '/Attempt1/anat/@default_subject/tess_cortex_mid_low_2000V.mat']);

% make brainstorm coordinate system consistent with MNI coordinates for
% plotting (in terms of axis directions)
cortex.Vertices = cortex.Vertices(:, [2 1 3]);
% in particular, the left-right axis needs to be flipped
cortex.Vertices(:, 1) = -cortex.Vertices(:, 1);

% same for highres_cortex
cortex_highres.Vertices = cortex_highres.Vertices(:, [2 1 3]);
cortex_highres.Vertices(:, 1) = -cortex_highres.Vertices(:, 1);

% same for highres_cortex
cortex_lowres.Vertices = cortex_lowres.Vertices(:, [2 1 3]);
cortex_lowres.Vertices(:, 1) = -cortex_lowres.Vertices(:, 1);


% calculate extrapolation from coarse to highres cortex, takes one minute
mi = []; in_normal_to_high = [];
for ii = 1:size(cortex_highres.Vertices, 1);
    %   ii
    [mi(ii) in_normal_to_high(ii)] = min(eucl(cortex_highres.Vertices(ii, :), cortex.Vertices));
end

mi = []; in_low_to_high = [];
for ii = 1:size(cortex_highres.Vertices, 1);
    %   ii
    [mi(ii) in_low_to_high(ii)] = min(eucl(cortex_highres.Vertices(ii, :), cortex_lowres.Vertices));
end

[xx, ia, ib] = intersect(cortex.Vertices, cortex_lowres.Vertices, 'rows');
[so ic] = sort(ib);
in_normal_to_low = ia(ic);

% load a 3-shell BEM model for the ICBM152 head with EGI Hydrocel 129 cap
% that I precomputed manually in Brainstorm (I can show you the necessary
% steps if needed)
headmodel = load([bs_folder '/Attempt1/data/Template/@default_study/headmodel_surf_openmeeg.mat']);

% make format compatible with my routines
leadfield = permute(reshape(headmodel.Gain, [], 3, nvox), [1 3 2]);

save([result_folder 'bs_results'], 'cortex', 'cortex_highres', 'cortex_lowres', 'leadfield', ...
    'in_normal_to_high', 'in_low_to_high', 'in_normal_to_low','-v7.3');

