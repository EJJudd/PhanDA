%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Phanerozoic GMST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/';
figname = 'Fig2_PhanerozoicGMST.png';
savefig = true;

% Select colormaps
% Main GMST colormap:
% Select colormap
cm = customcolormap(linspace(0,1,2),{'#515151','#FCFCFC'},75);
cm = [cm;0,0,0;flipud(cm)];

% Color for extinctions
plotcol = hex2rgb('#CA6702',1);
% Colormap for extinction uncertainty
gradcol = repmat(hex2rgb({'#FFFFFF','#EDCCAB','#FFFFFF'},1),2,1);
% Color for ice
icecol = [hex2rgb('#0a9396',1),.5];

% PART 1: LOAD DATA
% Directory details
assdate = '27Jul2023';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","Ndata","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")
ice = readtable("/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Data/SupplementalData/IceExtent_Macdonald2019.xlsx");
[icon,~,transperancy ] = imread("Extinction.png");

% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
sbCorr = [true, false];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SnowballCorr = ",string(sbCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
[GMST, dims] = combineruns(GMST,idx,1);
% Revise percentiles & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
GTS = combinestages(GTS,"GTS",Preferences,endsize,1);
% Calculate percentiles
prctiles = [5:1:95];
P = cell2mat(cellfun(@(x) prctile(x, prctiles), GMST, 'UniformOutput', false));
[prctiles,age_grid] = meshgrid(prctiles,GTS.Average);


% Calculate ice extent by stage
icestage = NaN(size(GTS,1),1);
for ii = 1:size(GTS,1)
    idx = find(ice.Age<=GTS.LowerBoundary(ii) & ice.Age>=GTS.UpperBoundary(ii));
    icestage(ii) = mean(ice.MaxLat(idx));
end
icestage(2:3) = 40;
icestage(isnan(icestage)) = 90;

%% PART 3: MAKE FIGURE
% (a) Initialize figure
fig=figure('Name','SummaryGMST','color','w'); 
fig.Units='inches';sPos = fig.Position;
fig.Position=[sPos(1),sPos(2),6.5*1.5,3.2*1.5];
fig.Units='pixels';    
pause(0.5)
% (b) Plot GMST
ax = axes('Position',[.05,.1125,.9,.875]); hold on
contourf(age_grid, P, prctiles, 152, 'LineColor', 'none')
colormap(cm)
ylim([6 49]);
ax.YTick = [5:5:40];
yl = ylim;
plot(GTS.Average,cell2mat(cellfun(@(x) median(x), GMST, 'UniformOutput', false)),'k-','LineWidth',2)
caxis([5 95])
% (c) Plot Extinctions
icon = double(icon);
icon(:,:,1) = plotcol(1);
icon(:,:,2) = plotcol(2);
icon(:,:,3) = plotcol(3);
timing = [444.5, GTS.LowerBoundary(GTS.Stage == "Famennian"), ...
    GTS.LowerBoundary(GTS.Stage == "Induan"), ...
    GTS.LowerBoundary(contains(GTS.Stage, "Hettangian")), ...
    GTS.UpperBoundary(GTS.Stage == "Maastrichtian")];
uncert = [1, 0, 0, 0, 0];
f = [1 2 5 4 1;2 3 6 5 2];
y = 7;
for ii = 1:numel(timing)
    e = imagesc([timing(ii)-8 timing(ii)+8], [y y+1.25], icon);
    set(e ,'AlphaData',transperancy);
    if uncert(ii) ~= 0
    v1 = [timing(ii)-uncert(ii),yl(1); ...
        timing(ii),yl(1); ...
        timing(ii)+uncert(ii),yl(1); ...
        timing(ii)-uncert(ii),yl(1)+.75; ...
        timing(ii),yl(1)+.75;...
        timing(ii)+uncert(ii),yl(1)+.75];
    v2 = [timing(ii)-uncert(ii),12; ...
        timing(ii),12; ...
        timing(ii)+uncert(ii),12; ...
        timing(ii)-uncert(ii),yl(2); ...
        timing(ii),yl(2);...
        timing(ii)+uncert(ii),yl(2)];
    patch('Faces',f,'Vertices',v1,...
      'FaceColor','interp','FaceVertexCData',gradcol,'edgecolor','none');
    patch('Faces',f,'Vertices',v2,...
      'FaceColor','interp','FaceVertexCData',gradcol,'edgecolor','none');
    end
    plot([timing(ii),timing(ii)],[yl(1) 6.75],'--','color',plotcol)
    plot([timing(ii),timing(ii)],[9 yl(2)],'--','color',plotcol)
end
% (d) Tidy axes
geologictimescale(0,GTS.LowerBoundary(size(P,1)),...
    'normal','reverse',gca,'standard','all','off',5.5,1)
yyaxis('left')
ax.FontSize = 11;ax.FontName = 'Arial';
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,...
    'FontWeight','bold')
% (e) Plot Ice extent
yyaxis right, ylim([-50,90]), ax.YTick = [40:10:90];
for ii = 1:size(icestage,1)
    x = GTS.UpperBoundary(ii);
    w = GTS.LowerBoundary(ii)-GTS.UpperBoundary(ii);
    y = icestage(ii);
    h = 90-icestage(ii);
    rectangle('Position',[x,y,w,h],'FaceColor',icecol,...
        'EdgeColor','none');
end
ax.YTickLabel = strsplit(num2str(ax.YTick),' ');
ax.YColor = icecol(1:3);
yl = ylabel(['Lat. ice extent (|',char(176),'|)'],'FontName','Arial','FontSize',13,...
    'FontWeight','bold','Color',icecol(1:3),'HorizontalAlignment','center');
yl.Position = [-13,65,-1];
% (f) add percentiles legend
pause(1)
leg=colorbar('north');
pause(1)
leg.Position=[.5865,.95,.225,.03];
leg.FontSize = 11;
pause(1)
xlabel(leg,'Percentile','FontName','Arial','FontSize',13,'FontWeight','bold')
% (g) Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end
