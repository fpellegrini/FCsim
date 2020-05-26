
load cm17
load('testdata.mat')
load('./processed_bs/bs_results.mat')
smooth_cortex = 0.35;
iroi_seed = 11; 
iroi_tar = 65; 

pos = cortex_highres.Vertices;

pos1 = pos(cortex_highres.Atlas(3).Scouts(iroi_seed).Vertices(1,1),:);
pos2 = pos(cortex_highres.Atlas(3).Scouts(iroi_tar).Vertices(1,1),:);

SurfSmoothIterations = ceil(300 * smooth_cortex * length(cortex_highres.Vertices) / 100000);
vc = tess_smooth(cortex_highres.Vertices, 1, SurfSmoothIterations, ...
    tess_vertconn(cortex_highres.Vertices, cortex_highres.Faces), 1);

data_in= gt; %zeros(1,length(cortex_highres.Curvature));
allplots_cortex_BS(cortex_highres1, data_in, [min(data_in) max(data_in)],...
    cm17a,'.', 0,['ground_thruth'],  ...
    {pos1, ...
    pos2});