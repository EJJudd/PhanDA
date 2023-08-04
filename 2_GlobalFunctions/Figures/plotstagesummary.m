function plotstagesummary(Yuse, Yeuse, proxytype, paleolat, paleolon,...
    taspost, gmst, stagelab, maplat, maplon, dem, tasprior, figdir)

load("GTS2020_PETM.mat")
load('HadCM3Coordinates.mat','Lat','Lon')
[Lon,Lat] = meshgrid(Lon,Lat);
LatWeight = cosd(Lat(:,1));

% Define additional requisite variables    
cm = flipud(customcolormap([linspace(0,.45,3),linspace(.55,1,5)],...
    {'#005F73','#0A9396','#94D2BD','#B4BE65','#E9D8A6','#FFB703','#CA6702','#9B2226'},50));
caxis_mean = [-45,45];
caxis_var = [0,7.5];
carbidx = contains(proxytype,'d18') & ~contains(proxytype,'p');

% Make figure
fig=figure('Name','StageSummary','NumberTitle','off','visible','on');
set(fig,'color','w');
fig.Units='inches';sPos = fig.Position;
fig.Position=[-15,sPos(2),8.5,11];
fig.Units='pixels';    
pause(0.5)

%Initialize tiles
t = tiledlayout(14,12,'Padding','none','TileSpacing','compact');
str = strsplit(stagelab,'_');
str = str2num(strrep(str(1),'S',''));
titlelab = sprintf('Stage %d: %s\n(%.2f - %.2f Ma)', ...
    str, GTS.Stage(str), GTS.UpperBoundary(str), GTS.LowerBoundary(str));
title(t,titlelab,'FontName','Arial','FontSize',15,'FontWeight','bold')
meantile1 = nexttile([4 6]);
meantile2 = nexttile([4 6]);
cmtile1 = nexttile([1 12]);
vartile1 = nexttile([4 6]);
vartile2 = nexttile([4 6]);
cmtile2 = nexttile([1 12]);
carbtile = nexttile([4 4]); hold on, box on
carbmaptile = nexttile([2 4]);
noncarbtile = nexttile(129,[4 4]); hold on, box on
noncarbmaptile = nexttile(149,[2 4]);

% Prior Mean
axes(meantile1)
ax = worldmap('World');
setm(ax,'meridianlabel','off','parallellabel','off')
pcolorm(Lat,Lon,median(tasprior,3))
caxis(meantile1,caxis_mean)
colormap(meantile1, cm)
gridm off, shading interp, framem('FlineWidth',1,'FEdgeColor','k') 
contourm(maplat.(stagelab),maplon.(stagelab),dem.(stagelab),...
    0,'k-','linewidth',1)
titlelab = sprintf('Median Prior \n (GMST = %.1f ± %.1f ^oC)', ...
    median(latweightgmst(tasprior)), ...
    std(latweightgmst(tasprior)));
title(meantile1,titlelab,'FontName','Arial','FontSize',12,'FontWeight','bold')

% Posterior Mean
axes(meantile2)
ax = worldmap('World');
setm(ax,'meridianlabel','off','parallellabel','off')
pcolorm(Lat,Lon,median(taspost,3))
caxis(meantile2,caxis_mean)
colormap(meantile2, cm)
gridm off, shading interp, framem('FlineWidth',1,'FEdgeColor','k') 
contourm(maplat.(stagelab),maplon.(stagelab),dem.(stagelab),...
    0,'k-','linewidth',1)
titlelab = sprintf('Median posterior \n (GMST = %.1f ± %.1f ^oC)', ...
    median(gmst), std(gmst));
title(meantile2,titlelab,'FontName','Arial','FontSize',12,'FontWeight','bold')

% Colorbar
c = colorbar(cmtile1,'south');
caxis(cmtile1,caxis_mean), colormap(cmtile1, cm)
c.FontSize = 10; c.FontName = 'Arial';
xlabel(c,'Air Temperature (^oC)','FontName','Arial','FontSize',12,'FontWeight','bold')
cmtile1.Visible = 'off';
pause(0.5)

% Prior Variance
axes(vartile1)
ax = worldmap('World');
setm(ax,'meridianlabel','off','parallellabel','off')
pcolorm(Lat,Lon,std(tasprior,[],3))
caxis(vartile1,caxis_var)
colormap(vartile1,cm)
gridm off, framem('FlineWidth',1,'FEdgeColor','k') 
contourm(maplat.(stagelab),maplon.(stagelab),dem.(stagelab),...
    0,'k-','linewidth',1)
title('Prior Std. Dev.','FontName','Arial','FontSize',12,'FontWeight','bold')

% Posterior Variance
axes(vartile2)
ax = worldmap('World');
setm(ax,'meridianlabel','off','parallellabel','off')
pcolorm(Lat,Lon,std(taspost,[],3))
caxis(vartile2,caxis_var)
colormap(vartile2,cm)
gridm off, framem('FlineWidth',1,'FEdgeColor','k') 
contourm(maplat.(stagelab),maplon.(stagelab),dem.(stagelab),...
    0,'k-','linewidth',1)
title('Posterior Std. Dev.','FontName','Arial','FontSize',12,'FontWeight','bold')

% Colorbar
axes(cmtile2)
c = colorbar(cmtile2,'south');
caxis(cmtile2,caxis_var)
colormap(cmtile2,cm)
c.FontSize = 10; c.FontName = 'Arial';
xlabel(c,'Air Temperature (^oC)','FontName','Arial','FontSize',12,'FontWeight','bold')
cmtile2.Visible = 'off';

% Carbonate plot & map
axes(carbmaptile)
ax = worldmap('World');
setm(ax,'meridianlabel','off','parallellabel','off')
gridm off, framem('FlineWidth',1,'FEdgeColor','k') 
contourm(maplat.(stagelab),maplon.(stagelab),dem.(stagelab),...
    0,'k-','linewidth',.5)
if sum(carbidx>0)
yyecrossplot(carbtile, carbmaptile, Yuse, Yeuse, proxytype, ...
    paleolat, paleolon, carbidx)
end

% Non-carbonate plot & map
axes(noncarbmaptile);
ax = worldmap('World');
setm(ax,'meridianlabel','off','parallellabel','off')
gridm off, framem('FlineWidth',1,'FEdgeColor','k') 
contourm(maplat.(stagelab),maplon.(stagelab),dem.(stagelab),...
    0,'k-','linewidth',.5)
if sum(~carbidx>0)
yyecrossplot(noncarbtile, noncarbmaptile, Yuse, Yeuse, ...
    proxytype, paleolat, paleolon, ~carbidx)
end

pause(0.5)
export_fig(fig,strcat(figdir,'/',stagelab,'.png'),'-p0.01','-m3','-nocrop','-painters')
close(fig)

end