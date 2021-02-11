function reg_param = fp_eloreta_crossval(signal_sensor,leadfield,nfold)

[n_sensors, n_voxels,n_dims] = size(leadfield);
%logarithmic spacing between 0.001 and 0.5
a = [0.001 0.5];
al = log10(a);
bl = linspace(al(1),al(2),10);
regs = 10.^bl;

signal_sensor1 = reshape(signal_sensor,n_sensors,[]);

for ifold = 1:nfold
        
    inds_test = randperm(n_sensors,round(n_sensors/nfold));
    inds_training = setdiff(1:n_sensors,inds_test);
    % inds = crossvalind('kfold',sensortrainingdata,4)
    count = 1;
    
    for ireg = regs

        A_eloreta = squeeze(mkfilt_eloreta_v2(leadfield(inds_training,:,:),ireg));
        
        A2 = reshape(A_eloreta, [numel(inds_training), n_dims * n_voxels]);
        
        signal_source = A2' * signal_sensor1(inds_training,:);
        
        leadfield1 = reshape(leadfield(inds_test,:,:), [numel(inds_test), n_dims*n_voxels]);
        signal_test = leadfield1*signal_source;
        signal_true = signal_sensor1(inds_test,:);
        
        mrsq(count,ifold) = mean(mean(signal_test - signal_true).^2)/mean(mean(signal_true.^2));
        
        count = count+1;
    end
end

reg_param = regs(find(mean(mrsq,2)==min(mean(mrsq,2))));
