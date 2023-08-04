%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GMST by Region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
savefig = true;

% Select colormap
cm = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);

% PART 1: LOAD DATA
% Directory details
assdate = '27Jul2023';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","LTG","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")

% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
sbCorr = [true, false];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",string(pHCorr))) & ...
    contains(ItName,strcat("SnowballCorr = ",string(sbCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);
LTG = combineruns(LTG,idx,2);

% Revise LTG, GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
LTG = combinestages(LTG,"LTG",Preferences,endsize,2);

% Define the regions
poles = abs(Lat)>=66.5;
tropics = abs(Lat)<=23.5;


%% FIGURE tropical comparison
% SAT and SST LTG
p = [5,16,50,84,95];
ltgsat = NaN(size(GMST,1),numel(p));
tropsat = NaN(size(GMST,1),numel(p));
meandim = 1;
for ii = 1:numel(GMST)
    xt = LTG{ii}; xt(~tropics,:) = NaN;
    xp = LTG{ii}; xp(~poles,:) = NaN;
    ltgsat(ii,:) = prctile(latweightgmst(xt,meandim)-latweightgmst(xp,meandim),p);
    tropsat(ii,:) = prctile(latweightgmst(xt,meandim),p);
end

fig = figure('Units','inches','Position',[-22,5,11.5,8.5],'Color','w');
tiledlayout(2,1,'Padding','none','TileSpacing','compact');
ax1 = nexttile; hold on, box on
fill([GTS.Average;flipud(GTS.Average)],[ltgsat(:,1);flipud(ltgsat(:,end))],'k','FaceColor','k','EdgeColor','none','FaceAlpha',.2);
fill([GTS.Average;flipud(GTS.Average)],[ltgsat(:,2);flipud(ltgsat(:,end-1))],'k','FaceColor','k','EdgeColor','none','FaceAlpha',.2);
plot(GTS.Average,ltgsat(:,3),'k-','LineWidth',2)
ylim([2 53]); ax1.YTick = [10:10:50];
geologictimescale(0,GTS.LowerBoundary(end),'normal','reverse',ax1,'standard','stages','off',9,2)
ylabel(['Latitudinal temperature gradient (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
text(ax1,482,50,'A','FontName','Arial','FontSize',15,'FontWeight','bold','color','k')


ax2 = nexttile; hold on, box on
fill([GTS.Average;flipud(GTS.Average)],[tropsat(:,1);flipud(tropsat(:,end))],'k','FaceColor',cm(end,:),'EdgeColor','none','FaceAlpha',.2);
fill([GTS.Average;flipud(GTS.Average)],[tropsat(:,2);flipud(tropsat(:,end-1))],'k','FaceColor',cm(end,:),'EdgeColor','none','FaceAlpha',.2);
plot(GTS.Average,tropsat(:,3),'-','LineWidth',2,'Color',cm(end,:))
ylabel(['Tropical temperature (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
geologictimescale(0,GTS.LowerBoundary(end),'normal','reverse',ax2,'standard','stages','off',9,2)
text(482,48,'B','FontName','Arial','FontSize',15,'FontWeight','bold','color','k')

box on
% Save figure
if savefig
    export_fig(gcf,[figdir,'/Supplemental/','SupFig_GMSTbyRegion.png'],'-p0.01','-m5')
end
