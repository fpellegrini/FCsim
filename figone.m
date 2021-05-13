function figone(height,width)
if ~exist('width','var')
    width = 8.9;
end
set(gcf,'Units','centimeters','PaperUnits','centimeters','Position',[5 5 width height],'PaperPosition',[5 5 width height],'PaperPositionMode','auto','PaperSize',[1+width 1+height])