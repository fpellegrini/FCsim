function fp_visualize_lf

patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

%%
for id = 7:length(patientID)
    
    load(sprintf('BF_Patient%s.mat',patientID{id}));
    clear L1 L
    L1 = inverse.MEG.L;
    
    for is=1:length(L1)
        L(is,:,:) = L1{is};
    end
    
    %% 1. One MEG sensor
    %     nchan(id) = size(L,2)
    isens = randi(size(L,2));
    clear pos
    pos = fp_getMNIpos(patientID{id});
    
%     dat = 1:size(L,1);
%     fp_plot_slices(dat,pos,id,isens,5,[1 size(L,1)])
    
    
        for idim = 1:3
            clear dat
            dat = squeeze(L(:,isens,idim)).*10^12;
            an(id,idim) =min(dat(:));
            ax(id,idim)=max(dat(:));
            fp_plot_slices(dat,pos,id,isens,idim,[100 (10^-5)+100])
                   % fp_data2nii(dat,pos,[],outname,id)
        end
end

%% 2. One voxel

id = 3;
ivox = 2665;

D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
D_ft = ftraw(D);
chanpos = D_ft.grad.chanpos;
chanpos(D.badchannels,:)=[];
chanpos(:,2) = -chanpos(:,2);
loc = mk_sensors_plane(chanpos(:,[2 1 3]));

figone(30,35)
for idim = 1:3
    clear dat
    dat = squeeze(L(ivox,:,idim)).*10^12;
%     dat(47)=1000
    
    subplot(2,2,idim)
    showfield_general(dat,loc);
    caxis([-2.5 2.5])
    
end
suptitle(['subject ' num2str(id), ', voxel ', num2str(ivox)])
outname = sprintf('sub_%d_voxel%d.png',id,ivox);
print(outname,'-dpng');
close all