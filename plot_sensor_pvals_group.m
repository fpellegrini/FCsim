

load('perm_sensor_beta_group.mat')

fileName = sprintf('redPLFP%s_off', '25');
D = spm_eeg_load(fileName);
 D_ft = ftraw(D);

pars.scale = [0 2.5];
loc = mk_sensors_plane(D_ft.grad.chanpos(1:125, [2 1 3]));
figure
figone(30,34)

%plot right channels



subplot(2,2,1)
showfield_general(-log10(p_l_abs), loc,pars); 
text(-0.5, 0.7,'left abs')

subplot(2,2,3)
showfield_general(-log10(p_l_im), loc,pars); 
text(-0.5, 0.7,'left imag')

subplot(2,2,2)
showfield_general(-log10(p_r_abs), loc,pars); 
text(-0.4, 0.7,'right abs')

subplot(2,2,4)
showfield_general(-log10(p_r_im), loc,pars); 
text(-0.5, 0.7,'right imag')
