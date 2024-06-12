%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% d18O-only DA Comparison %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
savefig = true;

% Select colormap
cm = hex2rgb({'#004F60';'#0a9396';'#ffb703';'#ca6702';'#9b2226'},1);

% PART 1: LOAD DATA
% Directory details
assdate = '27Jul2023';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","ItName")
GMST_all = GMST;
load([assdir,'/OutputWorkspaces/','Output_d18Oonly.mat'],"GMST","ItName","Ndata")
GMST_d18Oonly = GMST; clear GMST
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
GMST_all = combineruns(GMST_all,idx,1);
GMST_d18Oonly = combineruns(GMST_d18Oonly,idx,1);
Ndata = mean(Ndata,2);

% Revise GMST, & GTS to account for combined stages
endsize = size(GMST_all,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST_all = combinestages(GMST_all,"GMST",Preferences,endsize,1);
GMST_d18Oonly = combinestages(GMST_d18Oonly,"GMST",Preferences,endsize,1);
Ndata = combinestages(Ndata,"Ndata",Preferences,endsize,1);
Ndata = mean(Ndata,2);

paleozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Tremadocian"));
mesozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Selandian/Danian") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Induan"));
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));


%% MAKE FIGURE
figname = 'SupFig_d18OonlyComp.png';
fig = figure('Position',[-1248,491,479*2,420],'Color','w'); 
t = tiledlayout(1,5,'Padding','none','TileSpacing','compact');
ax1 = nexttile([1,3]); hold on
% Plot PhanDA GMST
GMST_all_median = cell2mat(cellfun(@(x) prctile(x,[5,16,50,84,95]), GMST_all, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[GMST_all_median(:,1);...
    flipud(GMST_all_median(:,end))],'k','FaceColor',[.85 .85 .85],'EdgeColor','none');
fill([GTS.Average;flipud(GTS.Average)],...
    [GMST_all_median(:,2);flipud(GMST_all_median(:,end-1))],'k','FaceColor',...
    [.65 .65 .65],'EdgeColor','none');
p1 = plot(GTS.Average,GMST_all_median(:,3),'k-','LineWidth',2);
% Plot d18O-only GMST
GMST_d18Oonly_median = cell2mat(cellfun(@(x) prctile(x,[5,16,50,84,95]), GMST_d18Oonly, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[GMST_d18Oonly_median(:,1);...
    flipud(GMST_d18Oonly_median(:,end))],'k','FaceColor',cm(1,:),'EdgeColor','none','FaceAlpha',.25);
fill([GTS.Average;flipud(GTS.Average)],...
    [GMST_d18Oonly_median(:,2);flipud(GMST_d18Oonly_median(:,end-1))],'k','FaceColor',...
    cm(1,:),'EdgeColor','none','facealpha',.5);
p1 = plot(GTS.Average,GMST_d18Oonly_median(:,3),'-','LineWidth',2,'color',cm(1,:));
ylim([0 45])
geologictimescale(0,GTS.LowerBoundary(end),...
    'normal','reverse',gca,'standard','all','off',9,2)
set(gca,'FontName','Arial','FontSize',11)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
text(450,8,"PhanDA",'FontWeight','bold','FontName','Arial','FontSize',11,'Color','k')
text(450,6.5,"\delta^{18}O-only DA",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(1,:))
text(480,43,"A",'FontWeight','bold','FontName','Arial','FontSize',15,'Color','k')

nexttile([1,2]), hold on, box on
% yneg = GMST_d18Oonly_median(:,3)-GMST_d18Oonly_median(:,1);
% ypos = GMST_d18Oonly_median(:,5)-GMST_d18Oonly_median(:,3);
% xneg = GMST_all_median(:,3)-GMST_all_median(:,1);
% xpos = GMST_all_median(:,5)-GMST_all_median(:,3);
% errorbar(GMST_all_median(:,3),GMST_d18Oonly_median(:,3),yneg,ypos,xneg,xpos,...
%     '.','Color',[.85,.85,.85],'LineWidth',.5,'CapSize',0)
yneg = GMST_d18Oonly_median(:,3)-GMST_d18Oonly_median(:,2);
ypos = GMST_d18Oonly_median(:,4)-GMST_d18Oonly_median(:,3);
xneg = GMST_all_median(:,3)-GMST_all_median(:,2);
xpos = GMST_all_median(:,4)-GMST_all_median(:,3);
errorbar(GMST_all_median(:,3),GMST_d18Oonly_median(:,3),yneg,ypos,xneg,xpos,...
    '.','Color',[.65,.65,.65],'LineWidth',1,'CapSize',0)
plot(GMST_all_median(paleozoic,3),GMST_d18Oonly_median(paleozoic,3),'ko','MarkerFaceColor',cm(4,:),'MarkerSize',10)
plot(GMST_all_median(mesozoic,3),GMST_d18Oonly_median(mesozoic,3),'ko','MarkerFaceColor',cm(3,:),'MarkerSize',10)
plot(GMST_all_median(Ndata(cenozoic)<=5,3),GMST_d18Oonly_median(Ndata(cenozoic)<=5,3),'wo','MarkerFaceColor',hex2rgb('#85C9CB',1),'MarkerSize',10)
plot(GMST_all_median(Ndata(cenozoic)>5,3),GMST_d18Oonly_median(Ndata(cenozoic)>5,3),'ko','MarkerFaceColor',cm(2,:),'MarkerSize',10)
plot(xlim,xlim,'k--','LineWidth',2)
text(15,36,"Cenozoic",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(2,:))
text(15,35,"Mesozoic",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(3,:))
text(15,34,"Paleozoic",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(4,:))
text(6,43.5,"B",'FontWeight','bold','FontName','Arial','FontSize',15,'Color','k')


% Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end

