%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Scenario Comparison  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sup. Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
savefig = true;        

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","Index","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")

% PART 2: PRE-TREAT DATA
% Select iterations to use
% (1a) pH: ens correction
pHCorr = "ens";
swCorr = ["snowball", "off"];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GpH_ens = combineruns(GMST,idx,1);
% (1b) pH: rec correction
pHCorr = "rec";
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GpH_rec = combineruns(GMST,idx,1);
% (1c) pH: no correction
pHCorr = "off";
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GpH_off = combineruns(GMST,idx,1);
% (2a) snowball: on correction
pHCorr = ["ens","rec"];
swCorr = "snowball";
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
Gsb_on = combineruns(GMST,idx,1);
% (2b) pH: no correction
swCorr = "off";
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
Gsb_off = combineruns(GMST,idx,1);
% (3a) R: low
swCorr = ["snowball", "off"];
rMeth = "low";
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
Gr_low = combineruns(GMST,idx,1);
% (3b) R: medium
rMeth = "medium";
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
Gr_med = combineruns(GMST,idx,1);
% (3b) R: high
rMeth = "high";
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
Gr_high = combineruns(GMST,idx,1);

% (4) All
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);



% Revise GMST & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GpH_ens = combinestages(GpH_ens,"GMST",Preferences,endsize,1);
GpH_rec = combinestages(GpH_rec,"GMST",Preferences,endsize,1);
GpH_off = combinestages(GpH_off,"GMST",Preferences,endsize,1);
Gsb_on = combinestages(Gsb_on,"GMST",Preferences,endsize,1);
Gsb_off = combinestages(Gsb_off,"GMST",Preferences,endsize,1);
Gr_low = combinestages(Gr_low,"GMST",Preferences,endsize,1);
Gr_med = combinestages(Gr_med,"GMST",Preferences,endsize,1);
Gr_high = combinestages(Gr_high,"GMST",Preferences,endsize,1);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);



%% Plot figure
figname = 'SupFig_Scenarios.png';
fig = figure('Position',[-1290,168,715,815],'Color','w');
t = tiledlayout(3,1,'Padding','none','TileSpacing','compact');
cm = hex2rgb({'#004F60';'#0a9396';'#ffb703';'#ca6702';'#9b2226'},1);
p = [16,50,84];
gmst = cell2mat(cellfun(@(x) prctile(x, p), GMST, 'UniformOutput', false));

% Panel 1: Snowball correction
ax = nexttile; hold on, box on
y1 = cell2mat(cellfun(@(x) prctile(x, p), Gsb_on, 'UniformOutput', false));
y2 = cell2mat(cellfun(@(x) prctile(x, p), Gsb_off, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[gmst(:,1);flipud(gmst(:,end))],'k',...
    'EdgeColor','none','FaceAlpha',.1)
fill([GTS.Average;flipud(GTS.Average)],[y1(:,1);flipud(y1(:,end))],'k',...
    'FaceColor',cm(1,:),'EdgeColor','none','FaceAlpha',.25)
fill([GTS.Average;flipud(GTS.Average)],[y2(:,1);flipud(y2(:,end))],'k',...
    'FaceColor',cm(end,:),'EdgeColor','none','FaceAlpha',.25)
plot(GTS.Average,y1(:,2),'LineWidth',2,'Color',cm(1,:))
plot(GTS.Average,y2(:,2),'LineWidth',2,'Color',cm(end,:))
plot(GTS.Average,gmst(:,2),'LineWidth',2,'Color','k')
text(10,12,'PhanDA Solution','FontName','Arial','FontWeight','bold',...
    'FontSize',11,'HorizontalAlignment','right')
text(10,10,'Including Snowball Earth δ^{18}O\fontsize{8}sw \fontsize{11} correction',...
    'FontName','Arial','FontWeight','bold','FontSize',11,'Color',cm(1,:),'HorizontalAlignment','right')
text(10,8,'Excluding Snowball Earth δ^{18}O\fontsize{8}sw \fontsize{11} correction',...
    'FontName','Arial','FontWeight','bold','FontSize',11,'Color',cm(end,:),'HorizontalAlignment','right')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylim([5.25 44.75])
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Tremadocian"),'normal',...
    'reverse',ax,'standard','stages','off',9,1)
text(480,42,'A','FontName','Arial','FontWeight','bold','FontSize',15)

% Panel 2: pH correction
ax = nexttile; hold on, box on
y1 = cell2mat(cellfun(@(x) prctile(x, p), GpH_off, 'UniformOutput', false));
y2 = cell2mat(cellfun(@(x) prctile(x, p), GpH_ens, 'UniformOutput', false));
y3 = cell2mat(cellfun(@(x) prctile(x, p), GpH_rec, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[gmst(:,1);flipud(gmst(:,end))],'k',...
    'EdgeColor','none','FaceAlpha',.1)
fill([GTS.Average;flipud(GTS.Average)],[y1(:,1);flipud(y1(:,end))],'k',...
    'FaceColor',cm(1,:),'EdgeColor','none','FaceAlpha',.25)
fill([GTS.Average;flipud(GTS.Average)],[y2(:,1);flipud(y2(:,end))],'k',...
    'FaceColor',cm(end-1,:),'EdgeColor','none','FaceAlpha',.25)
fill([GTS.Average;flipud(GTS.Average)],[y3(:,1);flipud(y3(:,end))],'k',...
    'FaceColor',cm(end,:),'EdgeColor','none','FaceAlpha',.25)
plot(GTS.Average,y1(:,2),'LineWidth',2,'Color',cm(1,:))
plot(GTS.Average,y2(:,2),'LineWidth',2,'Color',cm(end-1,:))
plot(GTS.Average,y3(:,2),'LineWidth',2,'Color',cm(end,:))
plot(GTS.Average,gmst(:,2),'LineWidth',2,'Color','k')
text(10,14,'PhanDA Solution','FontName','Arial','FontWeight','bold',...
    'FontSize',11,'HorizontalAlignment','right')
text(10,12,'No pH\fontsize{8}sw \fontsize{11} correction',...
    'FontName','Arial','FontWeight','bold','FontSize',11,'Color',cm(1,:),'HorizontalAlignment','right')
text(10,10,'Model prior pH\fontsize{8}sw \fontsize{11} correction',...
    'FontName','Arial','FontWeight','bold','FontSize',11,'Color',cm(end-1,:),'HorizontalAlignment','right')
text(10,8,'CO\fontsize{8}2 \fontsize{11}proxy pH\fontsize{8}sw \fontsize{11} correction',...
    'FontName','Arial','FontWeight','bold','FontSize',11,'Color',cm(end,:),'HorizontalAlignment','right')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylim([5.25 44.75])
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Tremadocian"),'normal',...
    'reverse',ax,'standard','stages','off',9,1)
box on
text(480,42,'B','FontName','Arial','FontWeight','bold','FontSize',15)

% Panel 3: R values
ax = nexttile; hold on, box on
y1 = cell2mat(cellfun(@(x) prctile(x, p), Gr_low, 'UniformOutput', false));
y2 = cell2mat(cellfun(@(x) prctile(x, p), Gr_med, 'UniformOutput', false));
y3 = cell2mat(cellfun(@(x) prctile(x, p), Gr_high, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[gmst(:,1);flipud(gmst(:,end))],'k',...
    'EdgeColor','none','FaceAlpha',.1)
fill([GTS.Average;flipud(GTS.Average)],[y1(:,1);flipud(y1(:,end))],'k',...
    'FaceColor',cm(1,:),'EdgeColor','none','FaceAlpha',.25)
fill([GTS.Average;flipud(GTS.Average)],[y2(:,1);flipud(y2(:,end))],'k',...
    'FaceColor',cm(end-1,:),'EdgeColor','none','FaceAlpha',.25)
fill([GTS.Average;flipud(GTS.Average)],[y3(:,1);flipud(y3(:,end))],'k',...
    'FaceColor',cm(end,:),'EdgeColor','none','FaceAlpha',.25)
plot(GTS.Average,y1(:,2),'LineWidth',2,'Color',cm(1,:))
plot(GTS.Average,y2(:,2),'LineWidth',2,'Color',cm(end-1,:))
plot(GTS.Average,y3(:,2),'LineWidth',2,'Color',cm(end,:))
plot(GTS.Average,gmst(:,2),'LineWidth',2,'Color','k')
text(10,14,'PhanDA Solution','FontName','Arial','FontWeight','bold',...
    'FontSize',11,'HorizontalAlignment','right')
text(10,12,'Low R',...
    'FontName','Arial','FontWeight','bold','FontSize',11,'Color',cm(1,:),'HorizontalAlignment','right')
text(10,10,'Medium R',...
    'FontName','Arial','FontWeight','bold','FontSize',11,'Color',cm(end-1,:),'HorizontalAlignment','right')
text(10,8,'High R',...
    'FontName','Arial','FontWeight','bold','FontSize',11,'Color',cm(end,:),'HorizontalAlignment','right')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylim([5.25 44.75])
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Tremadocian"),'normal',...
'reverse',ax,'standard','stages','off',9,1)
box on
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
text(480,42,'C','FontName','Arial','FontWeight','bold','FontSize',15)


if savefig
    export_fig(gcf,[figdir,'/Supplemental/',figname],'-p0.01','-m5')
end
