function plotltg(exprun,exps,expall,ltg,gmst,lab,threshold,cm,ax,printtext)

load("HadCM3Coordinates.mat","Lat")


for ii = 1:numel(exprun)
    color = cm(expall == exprun(ii),:);
    plot(ax,Lat,ltg(:,exps==exprun(ii)),'-','color',color)
end
plot(ax,Lat(:,1),mean(ltg,2),'k--','LineWidth',4)
ax.FontSize = 12;
xlim([-90 90])
ax.XTick = [-90:30:90];
if printtext
    text(0,mean(ylim)-.25*range(ylim),...
        sprintf('GMST = %.1f%sC\n%s\n%s',...
            median(gmst),char(176),string(splitlines(lab))),...
        'HorizontalAlignment','center','FontSize',8)
end

if ~threshold
    set(ax,'Color',[.5 .5 .5])
end

end