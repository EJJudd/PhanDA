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
PaleoClump.PaleoLatBand = string(PaleoClump.PaleoLatBand);
% Define lat bands
high = abs(Lat)>40 & abs(Lat)<=60;
mid = abs(Lat)>20 & abs(Lat)<=40;
low = abs(Lat)<=20;
PaleoClump.PaleoLatBand(abs(PaleoClump.PaleoLat)>40) = "high";
PaleoClump.PaleoLatBand(abs(PaleoClump.PaleoLat)>20 & abs(PaleoClump.PaleoLat)<=40) = "mid";
PaleoClump.PaleoLatBand(abs(PaleoClump.PaleoLat)<=20) = "low";

% Extract 
P = [5,16,50,84,95];
PhanDAsst = NaN(numel(PaleoClump.PublishedSST),numel(P));
Veizersst = NaN(numel(PaleoClump.PublishedSST),numel(P));
meandim = 1;
for ii = 1:numel(PaleoClump.PublishedSST)
    sidx = find(contains(GTS.Stage, PaleoClump.Stage(ii)));
    sltg = LTGsst{sidx};
    vltg = LTGsst_Veizer{sidx};
    if PaleoClump.PaleoLatBand(ii) == "low"
        sltg(~low,:) = NaN;
        vltg(~low,:) = NaN;
    elseif PaleoClump.PaleoLatBand(ii) == "mid"
        sltg(~mid,:) = NaN;
        vltg(~mid,:) = NaN;
    elseif PaleoClump.PaleoLatBand(ii) == "high"
        sltg(~high,:) = NaN;
        vltg(~high,:) = NaN;
    end
    PhanDAsst(ii,:) = prctile(latweightgmst(sltg,meandim),P);
    Veizersst(ii,:) = prctile(latweightgmst(vltg,meandim),P);
end

PaleoClump.Stage = string(PaleoClump.Stage);
PaleoClump.Stage(PaleoClump.Stage == "Aeronian") = "Aeronian/Rhuddanian";
PaleoClump.Stage(PaleoClump.Stage == "Rhuddanian") = "Aeronian/Rhuddanian";
stages = unique(PaleoClump.Stage);
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
lowestclump = NaN(numel(stages),2);
for ii = 1:numel(stages)
    %Find data within stage
    sidx = find(PaleoClump.Stage == stages(ii));
    % Identify unique lat bands
    lidx = string(unique(PaleoClump.PaleoLatBand(sidx)));
    % Index stage
    pidx = find(GTS.Stage == stages(ii));
    % Cycle through each lat band and plot the lowest clumped value
    for jj = 1:numel(lidx)
        slidx = find(PaleoClump.Stage == stages(ii) & PaleoClump.PaleoLatBand == lidx(jj));
        midx = find(PaleoClump.PublishedSST(slidx) == min(PaleoClump.PublishedSST(slidx)));
        v = [PhanDAsst(slidx(midx),3),PaleoClump.PublishedSST(slidx(midx))];
        u = [.5*range(PhanDAsst(slidx(midx),[2,4])),PaleoClump.PublishedUncertainty(slidx(midx))];
        if isnan(u(2))
            u(2) = 5;
        end
        lowestclump(ii,:) = v;
        h = plotellipses(v,u);
        h.FaceColor = [cm(ii,:),.25];
        if lidx == "low"
            marker = 'o';
        elseif lidx == "mid"
            marker = 'v';
        elseif lidx == "high"
            marker = 's';
        end
        plot(PhanDAsst(slidx(midx),3),PaleoClump.PublishedSST(slidx(midx)),...
            'MarkerFaceColor',cm(ii,:),'MarkerEdgeColor','none','marker',marker)
        h.EdgeColor = 'none';
    end
end
xlim([0 53])
ylim(xlim)
plot(xlim,ylim,'k--')
set(gca,'FontSize',11)
xlabel(['PhanDA SST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
ylabel(['Clumped SST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
text(1.5,50,'B','FontName','helvetica','FontSize',15,'FontWeight','bold')


% Panel 3
nexttile, hold on, box on
lowestclump = NaN(numel(stages),2);
for ii = 1:numel(stages)
    %Find data within stage
    sidx = find(PaleoClump.Stage == stages(ii));
    % Identify unique lat bands
    lidx = string(unique(PaleoClump.PaleoLatBand(sidx)));
    % Index stage
    pidx = find(GTS.Stage == stages(ii));
    % Cycle through each lat band and plot the lowest clumped value
    for jj = 1:numel(lidx)
        slidx = find(PaleoClump.Stage == stages(ii) & PaleoClump.PaleoLatBand == lidx(jj));
        midx = find(PaleoClump.PublishedSST(slidx) == min(PaleoClump.PublishedSST(slidx)));
        v = [Veizersst(slidx(midx),3),PaleoClump.PublishedSST(slidx(midx))];
        u = [.5*range(Veizersst(slidx(midx),[2,4])),PaleoClump.PublishedUncertainty(slidx(midx))];
        if isnan(u(2))
            u(2) = 5;
        end
        lowestclump(ii,:) = v;
        h = plotellipses(v,u);
        h.FaceColor = [cm(ii,:),.25];
        if lidx == "low"
            marker = 'o';
        elseif lidx == "mid"
            marker = 'v';
        elseif lidx == "high"
            marker = 's';
        end
        plot(Veizersst(slidx(midx),3),PaleoClump.PublishedSST(slidx(midx)),...
            'MarkerFaceColor',cm(ii,:),'MarkerEdgeColor','none','marker',marker)
        h.EdgeColor = 'none';
    end
end
xlim([0 53])
ylim(xlim)
plot(xlim,ylim,'k--')
set(gca,'FontSize',11)
xlabel(['Veizer & Prokoph correction SST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
ylabel(['Clumped SST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
text(1.5,50,'C','FontName','helvetica','FontSize',15,'FontWeight','bold')
plot(30,12,'ko')
text(32,12,['|latitude| ≤ 20',char(176)],'FontName','helvetica','FontSize',11)
plot(30,9,'kv')
text(32,9,['|latitude| = 21 - 40',char(176)],'FontName','helvetica','FontSize',11)
plot(30,6,'ks')
text(32,6,['|latitude| > 40',char(176)],'FontName','helvetica','FontSize',11)

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


