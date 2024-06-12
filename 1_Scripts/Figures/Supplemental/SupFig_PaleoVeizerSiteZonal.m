%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Sup Fig  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Veizer/Paleozoic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental/';
figname = 'SupFig_PaleoVeizer_Regional.png';
savefig = false;

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Documents/PhanDA/5_Outputs/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","LTGsst","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")
load("PhanerozoicCO2v9.mat", "PhanerozoicCO2")

% PART 2: PRE-TREAT DATA
% Select iterations to use
% Veizer run:
pHCorr = ["ens","rec"];
swCorr = "veizer";
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",swCorr)) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST_Veizer = combineruns(GMST,idx,1);
LTGsst_Veizer = combineruns(LTGsst,idx,2);
% Non-Veizer runs
swCorr = ["snowball","off"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",swCorr)) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);
LTGsst = combineruns(LTGsst,idx,2);

% Revise LTG, GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
LTGsst = combinestages(LTGsst,"LTG",Preferences,endsize,2);
GMST_Veizer = combinestages(GMST_Veizer,"GMST",Preferences,endsize,1);
LTGsst_Veizer = combinestages(LTGsst_Veizer,"LTG",Preferences,endsize,2);

% Load Paleozoic clumped data
PaleoClump = readtable('/Users/emilyjudd/Documents/PhanDA/3_GlobalFiles/ValidationData/PaleozoicClumped.xlsx');
PaleoClump.StandardDeviation(PaleoClump.DOI == "https://doi.org/10.1016/j.epsl.2018.02.001") = ...
    PaleoClump.StandardError(PaleoClump.DOI == "https://doi.org/10.1016/j.epsl.2018.02.001");
% Flag data with a temp > 60 or a std > 15
PaleoClump(PaleoClump.PublishedSST > 60, :) = [];
PaleoClump(PaleoClump.StandardDeviation > 15, :) = [];
% Average data by site/stage
PaleoClump.Stage = string(PaleoClump.Stage);
PaleoClump.LeadAuthor = string(PaleoClump.LeadAuthor);
PaleoClump.Stage(PaleoClump.Stage == "Aeronian") = "Aeronian/Rhuddanian";
PaleoClump.Stage(PaleoClump.Stage == "Rhuddanian") = "Aeronian/Rhuddanian";
% Convert to nearest grid cell
load("HadCM3Coordinates.mat", "Lat", "Lon")
PaleoClump.PaleoLat = Lat(dsearchn(Lat, PaleoClump.PaleoLat));
PaleoClump.PaleoLon = Lon(dsearchn(Lon, PaleoClump.PaleoLon));
sites = unique([PaleoClump.PaleoLat,PaleoClump.PaleoLon,string(PaleoClump.Stage)],'rows');
ClumpSites.PaleoLat = str2double(sites(:,1));
ClumpSites.PaleoLon = str2double(sites(:,2));
ClumpSites.Stage = sites(:,3);
ClumpSites = struct2table(ClumpSites);
ClumpSites.SSTmean = NaN(height(ClumpSites),1);
ClumpSites.SSTstd = NaN(height(ClumpSites),1);
N = 1000;
for ii = 1:height(ClumpSites)
    didx = find(PaleoClump.PaleoLat==ClumpSites.PaleoLat(ii) & ...
        PaleoClump.PaleoLon==ClumpSites.PaleoLon(ii) & ...
        PaleoClump.Stage==ClumpSites.Stage(ii));
    if all(~isnan(PaleoClump.StandardDeviation(didx)))
        d = NaN(numel(didx),N);
        for jj = 1:numel(didx)
            d(jj,:) = normrnd(PaleoClump.PublishedSST(didx(jj)),...
                PaleoClump.StandardDeviation(didx(jj)),1,N);
        end
        ClumpSites.SSTmean(ii) = mean(d,'all');
        ClumpSites.SSTstd(ii) = std(d,[],'all');
    elseif any(isnan(PaleoClump.StandardDeviation(didx))) & ...
            unique(PaleoClump.LeadAuthor(didx)) == "Came"
        ClumpSites.SSTmean(ii) = mean(PaleoClump.PublishedSST(didx));
        ClumpSites.SSTstd(ii) = std(PaleoClump.PublishedSST(didx));
    end
end

% Extract 
PhanDAsst = NaN(numel(ClumpSites.SSTmean),2);
Veizersst = NaN(numel(ClumpSites.SSTmean),2);
for ii = 1:numel(ClumpSites.Stage)
    lidx = find(abs(Lat-ClumpSites.PaleoLat(ii)) == min(abs(Lat-ClumpSites.PaleoLat(ii))));
    sidx = find(contains(GTS.Stage, ClumpSites.Stage(ii)));
    sltg = LTGsst{sidx};
    vltg = LTGsst_Veizer{sidx};
    PhanDAsst(ii,:) = [mean(sltg(lidx,:)),std(sltg(lidx,:))];
    Veizersst(ii,:) = [mean(vltg(lidx,:)),std(vltg(lidx,:))];
end
stages = unique(ClumpSites.Stage);
so = NaN(1,numel(stages));
for ii = 1:numel(stages)
    so(ii) = find(GTS.Stage == stages(ii));
end
[~,sidx] = sort(so);
stages = stages(sidx);


%% Make figure
fig = figure('Position',[396,267,876,555],'Color','w');
tiledlayout(2,3,'Padding','none','TileSpacing','compact');
cm = colormap(viridis(numel(stages)+1));
vc = cm(4,:);

% Panel 1
nexttile([1,3]), hold on
P = [5,16,50,84,95];
pgmst = cell2mat(cellfun(@(x) prctile(x, P), GMST, 'UniformOutput', false));
pgmst_veizer = cell2mat(cellfun(@(x) prctile(x, P), GMST_Veizer, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[pgmst(:,1);flipud(pgmst(:,5))],'k',...
    'FaceAlpha',.2,'FaceColor','k','EdgeColor','none');
fill([GTS.Average;flipud(GTS.Average)],[pgmst(:,2);flipud(pgmst(:,4))],'k',...
    'FaceAlpha',.2,'FaceColor','k','EdgeColor','none');
plot(GTS.Average,pgmst(:,3),'k-','LineWidth',1.5)
fill([GTS.Average;flipud(GTS.Average)],[pgmst_veizer(:,1);flipud(pgmst_veizer(:,5))],'k',...
    'FaceAlpha',.2,'FaceColor',vc,'EdgeColor','none');
fill([GTS.Average;flipud(GTS.Average)],[pgmst_veizer(:,2);flipud(pgmst_veizer(:,4))],'k',...
    'FaceAlpha',.2,'FaceColor',vc,'EdgeColor','none');
plot(GTS.Average,pgmst_veizer(:,3),'k-','LineWidth',1.5,'Color',vc)
ylim([-8 47])
set(gca,'YTick',[0:10:40])
geologictimescale(0,GTS.LowerBoundary(end),...
    'normal','reverse',gca,'standard','all','off',5.5,1,'helvetica',11,'k','k')
xlabel('Age (Ma)','FontName','helvetica','FontSize',13,'FontWeight','bold')
ylabel(['GMST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
text(480,43,'A','FontName','helvetica','FontSize',15,'FontWeight','bold')
text(225,3,'PhanDA','FontName','helvetica','FontSize',12,'FontWeight','bold')
text(225,0,'Veizer & Prokoph correction','FontName','helvetica','FontSize',12,'Color',vc,'FontWeight','bold')
plot(xlim,[14,14],'k--')
text(100,15.75,['Holocene GMST (14',char(176),'C)'],'FontName','helvetica','FontSize',11)

% Panel 2
nexttile, hold on, box on
for ii = 1:numel(stages)
    sidx = find(ClumpSites.Stage == stages(ii));
     for jj = 1:numel(sidx)
        v = [PhanDAsst(sidx(jj),1),ClumpSites.SSTmean(sidx(jj))];
        u = [PhanDAsst(sidx(jj),2),ClumpSites.SSTstd(sidx(jj))];
        h = plotellipses(v,u);
        h.FaceColor = [cm(ii,:),.25];
        h.EdgeColor = 'none';
    end
end
for ii = 1:numel(stages)
    sidx = find(ClumpSites.Stage == stages(ii));
    plot(PhanDAsst(sidx,1),ClumpSites.SSTmean(sidx),'ko',...
        'MarkerFaceColor',cm(ii,:),'MarkerEdgeColor','none')
end
xlim([0 60])
ylim(xlim)
plot(xlim,ylim,'k--')
set(gca,'FontSize',11)
xlabel(['PhanDA SST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
ylabel(['Clumped SST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
text(1.5,57,'B','FontName','helvetica','FontSize',15,'FontWeight','bold')
text(1.5,54,sprintf('r = %.2f',corr(ClumpSites.SSTmean,PhanDAsst(:,1))),'FontName','helvetica','FontSize',11,'FontWeight','bold')

% Panel 3
nexttile, hold on, box on
for ii = 1:numel(stages)
    sidx = find(ClumpSites.Stage == stages(ii));
     for jj = 1:numel(sidx)
        v = [Veizersst(sidx(jj),1),ClumpSites.SSTmean(sidx(jj))];
        u = [Veizersst(sidx(jj),2),ClumpSites.SSTstd(sidx(jj))];
        h = plotellipses(v,u);
        h.FaceColor = [cm(ii,:),.25];
        h.EdgeColor = 'none';
    end
end
for ii = 1:numel(stages)
    sidx = find(ClumpSites.Stage == stages(ii));
    plot(Veizersst(sidx,1),ClumpSites.SSTmean(sidx),'ko',...
        'MarkerFaceColor',cm(ii,:),'MarkerEdgeColor','none')
end
xlim([0 60])
ylim(xlim)
plot(xlim,ylim,'k--')
set(gca,'FontSize',11)
xlabel(['Veizer & Prokoph correction SST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
ylabel(['Clumped SST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
text(1.5,57,'C','FontName','helvetica','FontSize',15,'FontWeight','bold')
text(1.5,54,sprintf('r = %.2f',corr(ClumpSites.SSTmean,Veizersst(:,1))),'FontName','helvetica','FontSize',11,'FontWeight','bold')

% Panel 4
nexttile, hold on
for ii = 1:numel(stages)
    sidx = find(GTS.Stage==stages(ii));
    str = sprintf('{%s{%s}} (%.01f - %.01f Ma)','\bf',stages(ii),...
        GTS.LowerBoundary(sidx),GTS.UpperBoundary(sidx));
    text(.1,20-ii,str,'Color',cm(ii,:),...
        'FontName','helvetica','FontSize',12,'FontWeight','normal')
end
xlim([.075 1])
ylim([4 20])
set(gca,'Visible','off')


export_fig(gcf,[figdir,figname],'-p0.01','-m5')
