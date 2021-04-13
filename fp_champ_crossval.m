function reg_param = fp_champ_crossval(signal_sensor,leadfield,nfold)

tic

[n_sensors, n_voxels,n_dims] = size(leadfield);

regs = logspace(-2,0,15) * mean(diag(cov(signal_sensor(:,:)')));
regs(end) = []; 

signal_sensor1 = reshape(signal_sensor,n_sensors,[]);
L_perm = permute(leadfield,[1 3 2]);

for ifold = 1:nfold
        
    inds_test = randperm(n_sensors,round(n_sensors/nfold));
    inds_training = setdiff(1:n_sensors,inds_test);
    % inds = crossvalind('kfold',sensortrainingdata,4)
    count = 1;
    
    for ireg = regs        
        
        sigu = ireg*eye(n_sensors-numel(inds_test));        
        [~,~,w] = awsm_champ(signal_sensor1(inds_training, :), L_perm(inds_training, :), sigu, 200, 3, 2, 0);
        
        A_champ = real(reshape(w',size(leadfield(inds_training,:,:))));
        A2 = reshape(A_champ, [numel(inds_training), n_dims * n_voxels]);
        
        signal_source = A2' * signal_sensor1(inds_training,:);
        
        leadfield1 = reshape(leadfield(inds_test,:,:), [numel(inds_test), n_dims*n_voxels]);
        signal_test = leadfield1*signal_source;
        signal_true = signal_sensor1(inds_test,:);
        
        mrsq(count,ifold) = mean(mean(signal_test - signal_true).^2)/mean(mean(signal_true.^2));
        
        count = count+1;

    end
end

reg_param = regs(find(mean(mrsq,2)==min(mean(mrsq,2))));

toc
