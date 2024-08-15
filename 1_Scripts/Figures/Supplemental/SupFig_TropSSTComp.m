%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tropical SST  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sup. Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 1: LOAD DATA
% Directory details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
% PART 1: LOAD DATA
load([assdir,'/OutputWorkspaces/','Output.mat'],"ItName","LTGsst")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")

% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
swCorr = ["snowball", "off"];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",swCorr)) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
LTGsst = combineruns(LTGsst,idx,2);

% Revise LTG, GMST, & GTS to account for combined stages
endsize = size(LTGsst,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
LTGsst = combinestages(LTGsst,"LTGsst",Preferences,endsize,2);

% Define the regions
tropics = abs(Lat)<=30;

% Calculate mean tropic SST
tropsst = cell(size(LTGsst));
meandim = 1;
for ii = 1:numel(LTGsst)
    xt = LTGsst{ii}; xt(~tropics,:) = NaN;
    tropsst{ii} = (latweightgmst(xt,meandim))';
end
TropSST = cell2mat(cellfun(@(x) prctile(x,[5,16,50,84,95]), tropsst, 'UniformOutput', false));

% Load Grossman data
datadir = '/Users/emilyjudd/Documents/PhanDA/4_NonGlobalFiles/DataFiles';
load([datadir,'/GrossmanSST.mat'])

%% Make Figure
figname = 'SupFig_TropSSTComp.png';
fig = figure('Position',[-1248,491,479*2,420],'Color','w'); 
hold on, box on
% Plot PhanDA GMST
fill([GTS.Average;flipud(GTS.Average)],[TropSST(:,1);...
    flipud(TropSST(:,end))],'k','FaceColor',[.85 .85 .85],'EdgeColor','none');
fill([GTS.Average;flipud(GTS.Average)],...
    [TropSST(:,2);flipud(TropSST(:,end-1))],'k','FaceColor',...
    [.65 .65 .65],'EdgeColor','none');
p1 = plot(GTS.Average,TropSST(:,3),'k-','LineWidth',2);
plot(GrossmanSST.Age,GrossmanSST.SST,'-','color',hex2rgb('#CA6702',1),'LineWidth',2)
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Tremadocian"),'normal',...
    'reverse',gca,'standard','stages','off',8,1)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['Tropical SST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
text(483.5,19,'PhanDA','FontSize',13,'FontName','Arial')
text(483.5,17,'Grossman and Joachimski (2022)','FontSize',13,'FontName','Arial','Color',hex2rgb('#CA6702',1))

export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')