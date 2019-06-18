patientID = {'04'; '07'; '08'; '09'; '10';'11';'12';'18';'20';'22';'25'};

for id = 1
    clearvars -except patientID id
    D = spm_eeg_load(sprintf('redPLFP%s_off', patientID{id}));
    X = D(:,:,:);
    X1 = X;
   
    for ii = 1:size(X,1)
        clear o u 
        o = squeeze(X(ii,:,:));
        u = 10^(log10(range(o(:))));
        X1(ii,:,:)=X(ii,:,:)./u;
    end
      
    figure
    plot(X1(120,:,3)')
    hold on
    plot(X1(128,:,3)')
    
end