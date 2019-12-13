a =zeros(60,5);

for ii = [1 2 4:10 12:60]
    if size(varex{ii},1)<6
        ii
        a(ii,1:size(varex{ii},1)) = varex{ii};
    else
        a(ii,:)=varex{ii}(1:5);
    end
    plot(a')
    ylim([0.2 1])
    xlabel('Components')
    ylabel('Var explained (%)') 
    hold on 
end 
