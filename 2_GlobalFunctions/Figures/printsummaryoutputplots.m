function printsummaryoutputplots(OneDrive, InputDir, OutputDir, FigDir, ...
    OutputFilename, DArange, overwrite)

% Load paleogeographic maps
load([InputDir, '/YYe.mat'])
load([InputDir, '/Rvals.mat'])
load([OutputDir, '/', OutputFilename])
load([OutputDir, '/tasprior.mat'])

mapdirectory = [OneDrive '/Data/Paleogeographies/Scotese_PaleoDEMS'];
load('GTS2020_PETM.mat',"GTS")
load("ExpNo.mat", "ExpNo")
mapfilenames = selectmap(ExpNo.PlateAge, mapdirectory);
[~,maplat,maplon,dem] = scotesecontourlines(mapfilenames,GTS);
load('HadCM3Coordinates.mat','Lat','Lon')
[Lon,Lat] = meshgrid(Lon,Lat);

% Define additional requisite variables    
pHcorr = true;
snowballcorr = true;
Rmethod = 'percentile';
stagelab = strrep( strcat( ...
        'S', num2str((1:91)','%02.f'), '_', GTS.Stage(DArange)), ' ', '');
cm = flipud(customcolormap([linspace(0,.45,3),linspace(.55,1,5)],...
    {'#005F73','#0A9396','#94D2BD','#B4BE65','#E9D8A6','#FFB703','#CA6702','#9B2226'},50));
caxis_mean = [-40,40];
caxis_var = [0,25];



for a = DArange

    filename = [FigDir,'/ByStage/SummaryOutput_',char(stagelab(a)),'.png'];
    if overwrite || (~overwrite && ~exist(filename,'file'))

    % Extract Y Ye info
    [Yuse, Yeuse, Ruse, proxytype, paleolat, paleolon] = assembleYYeR( ...
        UPD.(stagelab{a}), Y.(stagelab{a}), Ye.(stagelab{a}), ...
        Rvals, Assumptions.(stagelab{a}), pHcorr, snowballcorr, Rmethod, a);
    carbidx = contains(proxytype,'d18') & ~contains(proxytype,'p');

    % Make figure
    fig=figure('Name','StageSummary','NumberTitle','off','visible','on');
    set(fig,'color','w');
    fig.Units='inches';sPos = fig.Position;
    fig.Position=[sPos(1),sPos(2),8.5,11];
    fig.Units='pixels';    
    pause(0.5)

    %Initialize tiles
    t = tiledlayout(14,12,'Padding','none','TileSpacing','compact');
    titlelab = sprintf('Stage %d: %s\n(%.2f - %.2f Ma)', ...
        a, GTS.Stage(a), GTS.UpperBoundary(a), GTS.LowerBoundary(a));
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
    pcolorm(Lat,Lon,taspriormean.(stagelab{a}))
    caxis(meantile1,caxis_mean)
    colormap(meantile1, cm)
    gridm off, shading interp, framem('FlineWidth',1,'FEdgeColor','k') 
    contourm(maplat.(stagelab{a}),maplon.(stagelab{a}),dem.(stagelab{a}),...
        0,'k-','linewidth',1)
    titlelab = sprintf('Prior Mean \n (GMST = %.1f ± %.1f ^oC)', ...
        gmstprior(a,1), gmstprior(a,2));
    title(meantile1,titlelab,'FontName','Arial','FontSize',12,'FontWeight','bold')

    % Posterior Mean
    axes(meantile2)
    ax = worldmap('World');
    setm(ax,'meridianlabel','off','parallellabel','off')
    pcolorm(Lat,Lon,TASmean.(stagelab{a}))
    caxis(meantile2,caxis_mean)
    colormap(meantile2, cm)
    gridm off, shading interp, framem('FlineWidth',1,'FEdgeColor','k') 
    contourm(maplat.(stagelab{a}),maplon.(stagelab{a}),dem.(stagelab{a}),...
        0,'k-','linewidth',1)
    titlelab = sprintf('Posterior Mean \n (GMST = %.1f ± %.1f ^oC)', ...
        mean(GMST(a,:)), std(GMST(a,:)));
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
    pcolorm(Lat,Lon,taspriorvar.(stagelab{a}))
    caxis(vartile1,caxis_var)
    colormap(vartile1,cm(.5*length(cm)+1:end,:))
    gridm off, framem('FlineWidth',1,'FEdgeColor','k') 
    contourm(maplat.(stagelab{a}),maplon.(stagelab{a}),dem.(stagelab{a}),...
        0,'k-','linewidth',1)
    title('Prior Variance','FontName','Arial','FontSize',12,'FontWeight','bold')

    % Posterior Variance
    axes(vartile2)
    ax = worldmap('World');
    setm(ax,'meridianlabel','off','parallellabel','off')
    pcolorm(Lat,Lon,TASvar.(stagelab{a}))
    caxis(vartile2,caxis_var)
    colormap(vartile2,cm(.5*length(cm)+1:end,:))
    gridm off, framem('FlineWidth',1,'FEdgeColor','k') 
    contourm(maplat.(stagelab{a}),maplon.(stagelab{a}),dem.(stagelab{a}),...
        0,'k-','linewidth',1)
    title('Posterior Variance','FontName','Arial','FontSize',12,'FontWeight','bold')

    % Colorbar
    axes(cmtile2)
    c = colorbar(cmtile2,'south');
    caxis(cmtile2,caxis_var)
    colormap(cmtile2,cm(.5*length(cm)+1:end,:))
    c.FontSize = 10; c.FontName = 'Arial';
    xlabel(c,'Air Temperature (^oC)','FontName','Arial','FontSize',12,'FontWeight','bold')
    cmtile2.Visible = 'off';

    % Carbonate plot & map
    axes(carbmaptile)
    ax = worldmap('World');
    setm(ax,'meridianlabel','off','parallellabel','off')
    gridm off, framem('FlineWidth',1,'FEdgeColor','k') 
    contourm(maplat.(stagelab{a}),maplon.(stagelab{a}),dem.(stagelab{a}),...
        0,'k-','linewidth',.5)
    yyecrossplot(carbtile, carbmaptile, Yuse, Yeuse, proxytype, ...
        paleolat, paleolon, carbidx)

    % Non-carbonate plot & map
    axes(noncarbmaptile);
    ax = worldmap('World');
    setm(ax,'meridianlabel','off','parallellabel','off')
    gridm off, framem('FlineWidth',1,'FEdgeColor','k') 
    contourm(maplat.(stagelab{a}),maplon.(stagelab{a}),dem.(stagelab{a}),...
        0,'k-','linewidth',.5)
    yyecrossplot(noncarbtile, noncarbmaptile, Yuse, Yeuse, ...
        proxytype, paleolat, paleolon, ~carbidx)

    pause(0.5)
    export_fig(filename,fig,'-p0.01','-m4','-nocrop','-painters')
    close(fig)

    end

end



end