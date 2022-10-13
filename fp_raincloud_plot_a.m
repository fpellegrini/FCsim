%% raincloud_plot - plots a combination of half-violin, boxplot,  and raw
% datapoints (1d scatter).
% Use as h = raincloud_plot(X, cl), where X is a data vector and cl is an
% RGB value. h is a cell array of handles for the various figure parts.
% Optional 3rd input argument 'density_type' can be 'ks' (default) or 'rash'.
% 'ks' uses matlab's 'ksdensity'. 'rash' uses 'rst_RASH' from Cyril
% Pernet's robust statistics toolbox (which must be on the matlab path).
% Based on https://micahallen.org/2018/03/15/introducing-raincloud-plots/
% Inspired by https://m.xkcd.com/1967/
% Written by Tom Marshall. www.tomrmarshall.com
% Thanks to Jacob Bellmund for some improvements
% Modified by Franziska Pellegrini


function [h, u] = fp_raincloud_plot_a(X, cl, box_on, bandwidth, density_type)

% if ~exist('density_type', 'var') | isempty(density_type)
%     density_type = 'ks'; % default is 'ks', can also be 'rash'
% end

if nargin < 3
    box_on = 1;
end

if nargin < 4
    bandwidth = [];
end

if nargin < 5
    density_type = 'ks';
end

%mean of the data 
mean_X = mean(X); 

%transform Y axis for better visibility 
base=100; X = -log(1-X*(1-1/base))/log(base);
mean_X_log = -log(1-mean_X*(1-1/base))/log(base);

% calculate kernel density
switch density_type
    case 'ks'
        [f, Xi, u] = ksdensity(X, 'bandwidth', bandwidth,'NumPoints',1000);
    case 'rash'
        try
            [Xi, f] = rst_RASH(X);
            u = NaN; % not sure how to handle this with RASH yet
        catch
            disp('you''ve specified density_type = ''RASH'', but something''s gone wrong.')
            disp('Have you downloaded Cyril Pernet''s robust stats toolbox?');
        end
end

% density plot
f_save = f;
f(Xi<0 | Xi>1) = [];
Xi(Xi<0 | Xi>1) = [];

h{1} = area(Xi, f); hold on
set(h{1}, 'FaceColor', cl);
set(h{1}, 'EdgeColor', 'none');
set(h{1}, 'LineWidth', 0.01);
plot([0 1],[0 0],'color',[1 1 1],'LineWidth',2)

% make some space under the density plot for the boxplot
yl = get(gca, 'YLim');
b = yl(2)/3;
set(gca, 'YLim', [-b*2 yl(2)]);

% width of boxplot
% wdth = b*0.5;
wdth = 0.5;

% jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;

% info for making boxplot
quartiles   = quantile(X, [0.25 0.75 0.5]);
iqr         = quartiles(2) - quartiles(1);
if iqr ~= 0 
    Xs          = sort(X);
    whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
    whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
    Y           = [quartiles whiskers];
else 
    Y = quartiles;
end

Y(end+1) = mean_X_log; 

% raindrops
h{2} = scatter(X, jit - 0.5);
%h{2} = scatter(X, jit - b/2);
h{2}.SizeData = 5;
h{2}.MarkerFaceColor = cl;
h{2}.MarkerEdgeColor = 'none';

if box_on
    % 'box' of 'boxplot'
    g = [0.3 0.3 0.3];
    h{3} = rectangle('Position', [Y(1) b/10-(wdth*0.5) Y(2)-Y(1) wdth]);   
    %h{3} = rectangle('Position', [Y(1) -b/2-(wdth*0.5) Y(2)-Y(1) wdth]);
    set(h{3}, 'EdgeColor', g)
    set(h{3}, 'LineWidth', 1);
    % could also set 'FaceColor' here as Micah does, but I prefer without
    
    % mean line    
    h{4} = line([Y(3) Y(3)], [b/10-(wdth*0.5) b/10+(wdth*0.5)], 'col', 'k', 'LineWidth', 1);
  h{5} = line([Y(end) Y(end)], [b/10-(wdth*0.5) b/10+(wdth*0.5)], 'col', 'r', 'LineWidth', 1);
    %h{4} = line([Y(3) Y(3)], [-b/2-(wdth*0.5) -b/2+(wdth*0.5)], 'col', 'r', 'LineWidth', 2);
    
    % whiskers
    if iqr ~= 0 
        h{6} = line([Y(2) Y(5)], [b/10 b/10], 'col',g, 'LineWidth', 1);
        h{7} = line([Y(1) Y(4)], [b/10 b/10], 'col',g, 'LineWidth', 1);
        %h{5} = line([Y(2) Y(5)], [-b/2 -b/2], 'col', 'k', 'LineWidth', 2);
        %h{6} = line([Y(1) Y(4)], [-b/2 -b/2], 'col', 'k', 'LineWidth', 2);
    end
end

xlim([0 1])
xticks1 = [0 0.5 0.7 0.9 0.95 0.97:0.01:1];
% xticks1 = [0 0.7 0.9 0.95 0.99 1];
xticks = -log(1-xticks1*(1-1/base))/log(base);
xTickLabels = xticks1;
% xTickLabels = -log(1-xTickLabels*(1-1/base))/log(base);
set(gca,'xtick',xticks,'xticklabels',xTickLabels);
