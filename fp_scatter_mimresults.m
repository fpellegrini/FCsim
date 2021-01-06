


load('./mim_sim1_auc.mat') %auc = 100 iterations x 18 pipelines x 3 measures 

figure
labs = {'mic','mim','mean icoh'};

for imim = 1:3 
    
    subplot(3,1,imim)
    
    data = auc(:,:,imim);

    for ii = 1:size(data,2)

        data1 = squeeze(data(:,ii));
        a = ii*3;

        ra = randi([-10 10],100,1)./20;
        scatter(ra+a,data1,'b')
        xlim([0 60])
        hold on
        plot([a-0.6 a+0.6],repmat(mean(data1),2,1),'r','linewidth',3)
        grid on
        hold on
    end
    
    xTicks = (1:size(data,2)).*3;
    xticklabels = {'1zs','2zs', '3zs','4zs','5zs',...
    '99% zs','90% zs','sumVox','baseline','99%corr','90%corr', ...
    '1 pc','2 pc', '3 pc', '4 pc', '5 pc', '99%', '90%'};
    ylabel(['auc ' labs{imim}])
    set(gca,'XTick', xTicks, 'XTickLabel',xticklabels)
    
end
