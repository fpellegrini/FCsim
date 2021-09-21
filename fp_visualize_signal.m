
ya = [-60 -40];
CS = tsdata_to_cpsd_fast(sensor_noise,fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

figone(12,10)
plot(u')
xlim([2 100])
ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)

xlabel('Frequency','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on 

saveas(gca,'~/Desktop/sensornoise.png','png')

close all
%%
ya = [-55 -15];

CS = tsdata_to_cpsd_fast(signal_sources',fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

figone(12,10)
plot(u')
xlim([2 100])
ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)

xlabel('Frequency','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on 

saveas(gca,'~/Desktop/interactive.png','png')

close all

%%
ya = [-55 -15];

CS = tsdata_to_cpsd_fast(noise_sources',fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

figone(12,10)
plot(u')
xlim([2 100])
ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)

xlabel('Frequency','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on 

saveas(gca,'~/Desktop/noninteractive.png','png')

close all


%%

ya = [-70 -49];
CS = tsdata_to_cpsd_fast(signal_sensor,fres,'WELCH');

for a=1:size(CS,1)
    u(a,:)=CS(a,a,:);
end

u=10*log10(u);

figone(12,10)
plot(u')
xlim([2 100])
ylim(ya)
xTicklabels = 0:10:50;
xTicks = 0:20:100;
set(gca,'XTick',xTicks,'XTickLabel',xTicklabels)

xlabel('Frequency','FontSize',16)
ylabel('Power (dB)','FontSize',16)
grid on 

saveas(gca,'~/Desktop/signal_sensor.png','png')

close all
