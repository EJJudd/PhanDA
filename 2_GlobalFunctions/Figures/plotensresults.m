function  plotensresults(gmst, ltg, idx, stagelab,...
    exprun, exps, explabs, exprunall, ...
    maplat, maplon, dem, tasprior, UPD, figdir)


cm1 = flipud(customcolormap(linspace(0,1,8),{'#59475C','#005f73','#0a9396','#b4be65','#ffb703','#ca6702','#bb3e03','#9b2226'},8));
cm2 = flipud(customcolormap([linspace(0,.45,3),linspace(.55,1,5)],...
    {'#005F73','#0A9396','#94D2BD','#B4BE65','#E9D8A6','#FFB703','#CA6702','#9B2226'},50));
load("GTS2020_PETM.mat")
load("ExpNo.mat", "ExpNo")
load("HadCM3Coordinates.mat")
LatWeight = cosd(Lat(:,1));
[Lon,Lat] = meshgrid(Lon,Lat);

fig = figure('Position',[-1872 1 1457 976],'Color','w'); 
tiledlayout(14,16,'Padding','none','TileSpacing','none');


% (1) Title
tlab = nexttile([2,6]);
text(.5,.75,strrep(stagelab,'_',' '),'FontSize',35,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle')
tlab.Visible = 'off';

% (2) Plot prior ltg
t = nexttile([6,4]); hold on, box on
for ii = 1:numel(exps{1})
    plot(Lat(:,1),squeeze(mean(tasprior(:,:,exprunall==exps{1}(ii)),2)),'-',...
        'color',cm1(ii,:))
end
yl = ylim;
y = linspace(yl(1)+.1*range(yl),yl(2)-range(yl)/2,numel(exps{1})+4);
for ii = 1:numel(exps{1})    
    text(0,y(ii),exps{1}(ii),'HorizontalAlignment','center','FontWeight','bold',...
        'Color',cm1(ii,:),'FontSize',12)
end
plot(Lat(:,1),squeeze(mean(median(tasprior,3),2)),'k--','LineWidth',4)
text(0,y(end-3),'mean of all','HorizontalAlignment','center','FontWeight','bold',...
        'Color','k','FontSize',12)
text(0,y(end),sprintf('%sGMST = %.1f^oC','\mu',median(latweightgmst(tasprior))),...
    'HorizontalAlignment','center','FontSize',15,'FontWeight','bold')
t.FontSize = 10; xlim([-90 90]), t.XTick = [-90:30:90];


% (3) Plot posterior ltg
t = nexttile([6,4]); hold on, box on
l = [];
g=[];
for ii = 1:size(idx,1)
    if idx{ii,3}
    for jj = 1:numel(exps{ii})
        color = cm1(exps{1} == exps{ii}(jj),:);
        plot(Lat(:,1),ltg(:,idx{ii,2}(exprun{ii}==exps{ii}(jj))),'-','color',color)
        l = [l,ltg(:,idx{ii,2}(exprun{ii}==exps{ii}(jj)))];
        g = [g;gmst(idx{ii,2}(exprun{ii}==exps{ii}(jj)))];
    end
    end
end
plot(Lat(:,1),median(l,2),'k--','LineWidth',4)
t.FontSize = 11;
xlim([-90 90])
t.XTick = [-90:30:90];
lab = ['All accepted solutions',newline,sprintf('(N = %.0f/35)',sum(cell2mat(idx(:,3))))];
text(0,mean(ylim)-.25*range(ylim),...
    sprintf('%s\n%s\n%sGMST = %.1f^o',...
        string(splitlines(lab)),'\mu',median(g)),...
    'HorizontalAlignment','center','FontSize',15,'FontWeight','bold')

    
% (4) Plot first assimilation
t = nexttile([2,2]); hold on, box on
plotltg(exprun{1},exps{1},exps{1},ltg(:,idx{1,2}),gmst(idx{1,2}),...
    explabs{1},idx{1,3},cm1,t)
t.FontSize = 10;

% (5) Maps
nexttile([4,6]);
ax = worldmap('World');setm(ax,'meridianlabel','off','parallellabel','off');
pcolorm(Lat,Lon,std(tasprior,[],3))
colormap(cm2), shading interp, framem k
caxis([0,7.5])
ca = caxis;
cb = colorbar(tlab,'Location','south');
caxis(tlab,ca)
ylabel(cb,'Prior STD_{SAT} (^oC)','FontSize',15,'FontWeight','bold')
contourm(maplat.(stagelab),maplon.(stagelab),dem.(stagelab),...
    0,'k-','linewidth',1)
fn = fieldnames(UPD.(stagelab));
for ii = 1:numel(fn)
    plotm(UPD.(stagelab).(fn{ii}).PaleoLat,UPD.(stagelab).(fn{ii}).PaleoLon,'kd','MarkerFaceColor','k')
end

% (6) Cycle through assimilations
for ii = 2:size(idx,1)
    t = nexttile([2,2]); hold on, box on
    plotltg(exprun{ii},exps{ii},exps{1},ltg(:,idx{ii,2}),...
        gmst(idx{ii,2}),explabs{ii},idx{ii,3},cm1,t)
    t.FontSize = 10;
end


export_fig(fig,strcat(figdir,'/',stagelab,'.png'),'-p0.01','-m3')
close(fig)

end