%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Sup Fig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Glacial LTG  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental/';
figname = 'SupFig_ColdhouseLTG.png';
savefig = true;

% Select colormap
cm = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);
cm_font = cm; cm_font(3,:) = [.5,.5,.5];

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Documents/PhanDA/5_Outputs/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","LTG","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")
load("PhanerozoicCO2v9.mat", "PhanerozoicCO2")

% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
swCorr = ["snowball", "off"];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",swCorr)) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);
LTG = combineruns(LTG,idx,2);
% Revise LTG, GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
LTG = combinestages(LTG,"LTG",Preferences,endsize,2);
mgmst = cell2mat(cellfun(@(x) nanmedian(x), GMST, 'UniformOutput', false));

%% Plot figure
cm = colormap(viridis(7));
cm = cm([1,3,4],:);
P = [16,50,84];
h = find(GTS.Stage=="Hirnantian");
    hltg = prctile(LTG{h},P,2);
    hltgfull = LTG{h};
    hice = sum(hltgfull<-10,2);
m = find(GTS.Stage=="Moscovian");
    mltg = prctile(LTG{m},P,2);
    mltgfull = LTG{m};
    mice = sum(mltgfull<-10,2);
t = find(GTS.Stage=="Tortonian");
    tltg = prctile(LTG{t},P,2);
    tltgfull = LTG{t};
    tice = sum(tltgfull<-10,2);
s = max([hice;mice;tice]);
fig = figure('Position',[939 402 410 420],'Color','w'); 
tiledlayout(7,1,'Padding','none','TileSpacing','compact');
ax1 = nexttile;  hold on, box on
y = 0;
w = 2.5;
hi = 1;
for ii = numel(Lat):-1:floor(numel(Lat)/2)
    x = Lat(ii);
    rectangle(ax1,'Position',[x,y,w,hi],'FaceColor',[cm(1,:),hice(ii)/s],...
        'EdgeColor','none')
    rectangle(ax1,'Position',[x,y-1,w,hi],'FaceColor',[cm(2,:),mice(ii)/s],...
        'EdgeColor','none')
    rectangle(ax1,'Position',[x,y-2,w,hi],'FaceColor',[cm(3,:),tice(ii)/s],...
        'EdgeColor','none')
end
for ii = 1:floor(numel(Lat)/2)
    x = Lat(ii+1);
    rectangle(ax1,'Position',[x,y,w,hi],'FaceColor',[cm(1,:),hice(ii)/s],...
        'EdgeColor','none')
    rectangle(ax1,'Position',[x,y-1,w,hi],'FaceColor',[cm(2,:),mice(ii)/s],...
        'EdgeColor','none')
    rectangle(ax1,'Position',[x,y-2,w,hi],'FaceColor',[cm(3,:),tice(ii)/s],...
        'EdgeColor','none')
end
set(ax1,'XLim',[-90,90],'XTick',[-90:30:90],'YTick',[])
text(ax1,0,-.5,'latitudinal ice extent','VerticalAlignment','middle',...
    'HorizontalAlignment','center','FontName','helvetica','FontSize',12,...
    'FontWeight','bold')
ax2 = nexttile([6,1]);  hold on, box on
fill([Lat;flipud(Lat)],[hltg(:,1);flipud(hltg(:,3))],'k','EdgeColor','none',...
    'FaceColor',cm(1,:),'FaceAlpha',.25)
fill([Lat;flipud(Lat)],[mltg(:,1);flipud(mltg(:,3))],'k','EdgeColor','none',...
    'FaceColor',cm(2,:),'FaceAlpha',.25)
fill([Lat;flipud(Lat)],[tltg(:,1);flipud(tltg(:,3))],'k','EdgeColor','none',...
    'FaceColor',cm(3,:),'FaceAlpha',.25)
plot(Lat,hltg(:,2),'-','Color',cm(1,:),'LineWidth',1.5)
plot(Lat,mltg(:,2),'-','Color',cm(2,:),'LineWidth',1.5)
plot(Lat,tltg(:,2),'-','Color',cm(3,:),'LineWidth',1.5)
set(gca,'XLim',[-90,90],'XTick',[-90:30:90])
str = sprintf('{%s{%s}}%s(%.01f - %.01f Ma)%sGMST = %.1f%sC','\bf',GTS.Stage(h),...
        newline,GTS.LowerBoundary(h),GTS.UpperBoundary(h),newline,mgmst(h),char(176));
text(0,5,str,'HorizontalAlignment','center','FontName','helvetica',...
    'FontSize',12,'Color',cm(1,:))
str = sprintf('{%s{%s}}%s(%.01f - %.01f Ma)%sGMST = %.1f%sC','\bf',GTS.Stage(m),...
        newline,GTS.LowerBoundary(m),GTS.UpperBoundary(m),newline,mgmst(m),char(176));
text(0,-7.5,str,'HorizontalAlignment','center','FontName','helvetica',...
    'FontSize',12,'Color',cm(2,:))
str = sprintf('{%s{%s}}%s(%.01f - %.01f Ma)%sGMST = %.1f%sC','\bf',GTS.Stage(t),...
        newline,GTS.LowerBoundary(t),GTS.UpperBoundary(t),newline,mgmst(t),char(176));
text(0,-20,str,'HorizontalAlignment','center','FontName','helvetica',...
    'FontSize',12,'Color',cm(3,:))
ylabel(['Temperature (',char(176),'C)'],'FontName','Arial','FontSize',15,'FontWeight','bold')
xlabel(['Latitude (',char(176),')'],'FontName','Arial','FontSize',15,'FontWeight','bold')

export_fig(gcf,[figdir,figname],'-p0.01','-m5')


