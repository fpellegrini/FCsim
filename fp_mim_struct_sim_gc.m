function fp_mim_struct_sim_gc(params)

DIROUT = '/home/bbci/data/haufe/Franziska/data/gc_sim/';
if ~exist(DIROUT);mkdir(DIROUT); end
DIROUT1 = '/home/bbci/data/haufe/Franziska/data/gc_save/';
if ~exist(DIROUT1);mkdir(DIROUT1); end

if params.ip==7 || params.ip==8
    params_save = params; 
    load(sprintf('%s/gc_CS/%d.mat',DIROUT1,params.iit));
    params = params_save; 
    clear params_save
else
    
    if params.ip==5 || params.ip ==4
        params_save = params;
        load(sprintf('%s/gc_sig/%d.mat',DIROUT1,params.iit));
        params= params_save;
        clear params_save 
    else
        
        %% signal generation
        tic
        % ROI labels
        % In D.sub_ind_roi, there are the randomly
        % selected voxels of each region
        fprintf('Getting atlas positions... \n')
        D = fp_get_Desikan(params.iReg);
        
        %signal generation
        fprintf('Signal generation... \n')
        [sig,brain_noise,sensor_noise,gt,L,iroi_seed, iroi_tar,D, fres, n_trials] = ...
            fp_generate_mim_signal_gc(params,D,DIROUT1);
     
        if params.ip==1
            dir1 =  sprintf('%s/gc_sig/',DIROUT1);
            if ~exist(dir1); mkdir(dir1); end
            outname = sprintf('%s/gc_sig/%d.mat',DIROUT1,params.iit);
            save(outname,'-v7.3')
        end
        
        t.signal = toc;
    end
       
    tic
    %combine noise sources
    noise = params.iss*brain_noise + (1-params.iss)*sensor_noise;
    noise = noise ./ norm(noise, 'fro');
    %combine signal and noise
    signal_sensor1 = params.isnr*sig + (1-params.isnr)*noise;
    signal_sensor = signal_sensor1 ./ norm(signal_sensor1, 'fro');
    
    %reshape
    signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
    
    t.snr = toc;
    
    if params.ip==1
       dir1 =  sprintf('%s/gc_CS/',DIROUT1);
       if ~exist(dir1); mkdir(dir1); end
       outname = sprintf('%s/gc_CS/%d.mat',DIROUT1,params.iit);
       save(outname,'-v7.3')
    end
end

%% leadfield
tic
L3 = L(:, D.ind_cortex, :);
L_backward = L3; 
ni = size(L_backward,3);

%construct source filter
if strcmp(params.ifilt,'e')
    reg_param = fp_eloreta_crossval(signal_sensor,L_backward,5);
    A = squeeze(mkfilt_eloreta_v2(L_backward,reg_param));
    A = permute(A,[1, 3, 2]);
    
    
elseif strcmp(params.ifilt,'l')
    cov_ = cov(data(:,:)');
    reg = 0.05*trace(cov_)/length(cov_);
    Cr = cov_ + reg*eye(size(cov_,1));
    
    [~, A] = lcmv(Cr, L_backward, struct('alpha', 0, 'onedim', 0));
    A = permute(A,[1, 3, 2]);
end

t.filter = toc;

%% calculate GC

if params.ip ==1
    
    %pca pipeline ('all' 8 pipelines + baseline)
    zs=1;
    [TRGC, GC, DIFFGC,to_save,t] = fp_get_gc(A,signal_sensor, D,'all',zs,t);
    fprintf('pipelines calculated')
    
    %% without ZS standardisation
    fprintf('Calculating fix, max and percent without ZS standardisation... \n')
    zs=0;
    for ii = 1:5
        [TRGC_fixed_zs0{ii}, GC_fixed_zs0{ii}, DIFFGC_fixed_zs0{ii}, to_save_fixed_zs0{ii},t1] = ...
            fp_get_gc(A,signal_sensor, D,ii,zs,[]);
        t.zs0fixed{ii} = t1.gc;
    end
    [TRGC_max_zs0, GC_max_zs0,  DIFFGC_max_zs0, to_save_max_zs0,t1] = ...
        fp_get_gc(A,signal_sensor, D,'max',zs,[]);
    t.zs0max = t1.gc;
    [TRGC_percent_zs0, GC_percent_zs0, DIFFGC_percent_zs0, to_save_percent_zs0,t1] = ...
        fp_get_gc(A,signal_sensor, D,'percent',zs,[]);
    t.zs0percent = t1.gc;

    TRGC.fixed_zs0 = TRGC_fixed_zs0;
    GC.fixed_zs0 = GC_fixed_zs0; 
    to_save.fixed_zs0 = to_save_fixed_zs0;
    DIFFGC.fixed_zs0 = DIFFGC_fixed_zs0; 

    TRGC.max_zs0 = TRGC_max_zs0; 
    GC.max_zs0 = GC_max_zs0; 
    to_save.max_zs0 = to_save_max_zs0;
    DIFFGC.max_zs0 = DIFFGC_max_zs0; 

    TRGC.percent_zs0 = TRGC_percent_zs0; 
    GC.percent_zs0 = GC_percent_zs0; 
    to_save.percent_zs0 = to_save_percent_zs0;
    DIFFGC.percent_zs0 = DIFFGC_percent_zs0;

    clear TRGC_max_zs0 GC_max_zs0 to_save_max_zs0 DIFFGC_max_zs0 TRGC_percent_zs0 ...
        GC_percent_zs0 to_save_percent_zs0 DIFFGC_percent_zs0 TRGC_fixed_zs0 ...
        GC_fixed_zs0 to_save_fixed_zs0 DIFFGC_fixed_zs0
else
    %reduced pca pipelines: zs0 of all fixed pips and of 90 and 99 %,
    %baseline 
    zs=0;
    [TRGC,GC,DIFFGC,to_save,t] = fp_get_gc_reduced(A,signal_sensor, D,'all',zs,t);
    fprintf('pipelines calculated')
end

to_save.t = t;

%% save performance and baseline
fprintf('Saving... \n')
tic
outname = sprintf('%sgc_%s.mat',DIROUT,params.logname);
save(outname,'-v7.3')
t.saving = toc;

