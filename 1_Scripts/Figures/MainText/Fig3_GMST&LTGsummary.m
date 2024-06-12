%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Climate states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
figname = 'Fig3_GMST&LTGsummary.png';
savefig = false;

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
%swCorr = ["snowball", "veizer","off"];
swCorr = "veizer";
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

% Define the regions
poles = abs(Lat)>=66.5;
tropics = abs(Lat)<=23.5;
% subdivide GMST by percentile
prctiles = [20:20:100];
mgmst = cell2mat(cellfun(@(x) nanmedian(x), GMST, 'UniformOutput', false));
P = prctile(mgmst, prctiles);
CSidx.ih = find(mgmst<=P(1));
CSidx.ch = find(mgmst>P(1) & mgmst<=P(2));
CSidx.tr = find(mgmst>P(2) & mgmst<=P(3));
CSidx.gh = find(mgmst>P(3) & mgmst<=P(4));
CSidx.hh = find(mgmst>P(4));
rmgmst = mgmst+abs(min(mgmst))+1e-5;
% subdivide by era
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
mesozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Selandian/Danian") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Induan"));
paleozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Tremadocian"));
% proportion of time spent in each state in each era
fn = fieldnames(CSidx);
for ii = 1:numel(fn)
	CStime.(fn{ii}) = [sum(GTS.LowerBoundary(CSidx.(fn{ii})) - GTS.UpperBoundary(CSidx.(fn{ii}))); ...
        sum(GTS.LowerBoundary(intersect(CSidx.(fn{ii}),cenozoic)) - GTS.UpperBoundary(intersect(CSidx.(fn{ii}),cenozoic))); ...
        sum(GTS.LowerBoundary(intersect(CSidx.(fn{ii}),mesozoic)) - GTS.UpperBoundary(intersect(CSidx.(fn{ii}),mesozoic))); ...
        sum(GTS.LowerBoundary(intersect(CSidx.(fn{ii}),paleozoic)) - GTS.UpperBoundary(intersect(CSidx.(fn{ii}),paleozoic)))];
end
CStime.Total = [CStime.ih,CStime.ch,CStime.tr,CStime.gh,CStime.hh];

%% PART 3: MAKE FIGURE
% (a) Initialize figure
fig = figure('Units','inches','Position',[2,2,6.5*1.5,4.3*1.5],'Color','w');

% (b) Plot color stripes
ax1 = subplot(3,6,[1:6]); hold on
ax1.Position = [.025, ax1.Position(2), .95, ax1.Position(4)+.05];
DArange = [1:size(GMST,1)];
% y = 7;
% h = 40;
% for ii = 1:numel(fn)
%     for jj = 1:numel(CSidx.(fn{ii}))
%         x = GTS.UpperBoundary(CSidx.(fn{ii})(jj));
%         w = GTS.LowerBoundary(CSidx.(fn{ii})(jj))-GTS.UpperBoundary(CSidx.(fn{ii})(jj));
%         rectangle('Position',[x,y,w,h],'EdgeColor','none','FaceColor',cm(ii,:))
%     end
% end
% ylim([8,42])
y = 0;
r = max(mgmst);
hold on
for ii = 1:numel(fn)
    for jj = 1:numel(CSidx.(fn{ii}))
        x = GTS.UpperBoundary(CSidx.(fn{ii})(jj));
        w = GTS.LowerBoundary(CSidx.(fn{ii})(jj))-GTS.UpperBoundary(CSidx.(fn{ii})(jj));
        h = abs(mgmst(CSidx.(fn{ii})(jj)))/r;
        rectangle('Position',[x,y,w,h],'EdgeColor','none','FaceColor',cm(ii,:))
    end
end
ylim([0 1.05])
geologictimescale(0,GTS.LowerBoundary(DArange(end)),...
    'normal','reverse',gca,'standard','stages','off',4,2)  
ax1 = gca; ax1.FontSize = 11; ax1.FontName = 'helvetica'; ax1.YTick = [];
ylabel('Climate state','FontName','helvetica','FontSize',13,'FontWeight','bold')
xlabel('Age (Ma)','FontName','helvetica','FontSize',13,'FontWeight','bold')
% text(482,38,'A','FontName','helvetica','FontSize',15,'FontWeight','bold','color','w')
text(ax1,484,.975,'A','FontName','helvetica','FontSize',15,'FontWeight','bold','color','k')
ax1.YRuler.Axle.Visible = 'off';

% (c) Distributions
ax2 = subplot(3,2,[3,5]);
ax2.Position = [.06,.075,.4,.4];
hold on
xedge = [floor(min(mgmst))-8:.25:ceil(max(mgmst))+8]';
y = ksdensity(mgmst,xedge);
xtext = [12,16,24.5,32.75,37];
ytext = [.0375,.05375,.07,.05375,.0375]';
statelab = ["Coldhouse","Coolhouse","Transitional","Warmhouse","Hothouse"];
for ii = 1:numel(fn)
    if ii == 1
        idx = xedge <= P(ii);
        temptext = sprintf("(%.0f-%.0f%sC)",round(min(mgmst)),round(P(ii)),char(176));
    elseif ii == numel(fn)
        idx = xedge > P(ii-1);
        temptext = sprintf("(%.0f-%.0f%sC)",round(P(ii-1)),round(max(mgmst)),char(176));
    else
        idx = xedge > P(ii-1) & xedge <= P(ii);
        temptext = sprintf("(%.0f-%.0f%sC)",round(P(ii-1)),round(P(ii)),char(176));
    end
    plotx = [xedge(idx);flipud(xedge(idx))];
    ploty = [y(idx);zeros(sum(idx),1)];
    fill(plotx,ploty,cm(ii,:),'FaceAlpha',1,'EdgeColor','none');
   text(xtext(ii),ytext(ii),statelab(ii),'FontWeight','bold','FontName',...
       'helvetica','FontSize',11,'HorizontalAlignment','center','Color',cm_font(ii,:))
   text(xtext(ii),ytext(ii)-.005,temptext,'FontName','helvetica','FontSize',...
       11,'HorizontalAlignment','center','Color',cm_font(ii,:))
end
ylim([0 .072]),xlim([3 44])
ax2.FontName = 'helvetica'; ax2.FontSize = 11;
xlabel(['GMST (',char(176),'C)'],'FontName','helvetica','FontSize',13,...
    'FontWeight','bold')
ylabel('Density','FontName','helvetica','FontSize',13,'FontWeight','bold')
text(3.75,.0695,'B','FontName','helvetica','FontSize',15,'FontWeight','bold')
ax2.YTick = [0:.01:.06];

% (d) Add panel C text and axes
ax3a = axes('Position',[.6,.1,.1,.1])
ax3a.Position = [.09,.395,.225,.245];
text(ax3a,1.05,.65, sprintf('Proportion of time\nspent in each state'),...
    'FontName','helvetica','FontSize',11,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','middle')
text(ax3a,1,.85,'C','FontName','helvetica','FontSize',15,'FontWeight','bold')
ax3a.Visible = 'off';
% (e) Pie chart
ax3b = axes('Position',[.6,.1,.1,.1])
ax3b.Position = [.28,.395,.245,.245];
h = pie(CStime.Total(1,:));
for ii = 2:2:2*numel(CStime.Total(1,:))
    h(ii).Position = 0.4*h(ii).Position;
    if ii ~= 6
        h(ii).Color = 'w';
    end
    h(ii).FontName = 'helvetica';
    h(ii).FontWeight = 'bold';
    h(ii).FontSize = 11;
    h(ii-1).EdgeColor = 'w';
end
colormap(cm)
% manually adjust labels to be in center of each slice
h(2).Position = [-.05,.475,0];
h(4).Position = [-.375,.115,0];
h(6).Position = [-.125,-.35,0];
h(8).Position = [.375,-.15,0];
h(10).Position = [.1,.35,0];

% (d) Plot latitudinal gradients
ax4 = subplot(3,6,[10:12,16:18]);
ax4.Position = [.55,.075,.425,.54];
hold on, box on
for ii = 1:numel(fn)
    % 1 sigma bound
    y=[];
    for jj = 1:numel(CSidx.(fn{ii}))
        y = [y, LTG{CSidx.(fn{ii})(jj)}];
    end
    f = fill([Lat;flipud(Lat)],[prctile(y,16,2);...
        flipud(prctile(y,84,2))],cm(ii,:),...
        'FaceAlpha',.25,'EdgeColor','none');
    plot(Lat,median(y,2),'-','color',cm_font(ii,:),'LineWidth',3)
end
% Tidy plot
xlim([-90,90]);
ax4.XTick = [-90:30:90];
ax4.FontName = 'helvetica';
ax4.FontSize = 11;
xlabel('Latitude','FontName','helvetica','FontSize',13,'FontWeight','bold')
ylabel(['Air temperature (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
ylim([-50 46])
text(-86,41.25,'D','FontName','helvetica','FontSize',15,'FontWeight','bold','color','k')
ax4.Color = 'none';

% (e) Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end


%% ALT
% (a) Initialize figure
fig = figure('Units','inches','Position',[2,2,6.5*1.5,4.6*1.5],'Color','w');

% (b) Plot color stripes
ax1 = subplot(3,6,[1:6]); hold on
ax1.Position = [.025, ax1.Position(2), .95, ax1.Position(4)+.05];
DArange = [1:size(GMST,1)];
y = 7;
h = 40;
hold on
for ii = 1:numel(fn)
    for jj = 1:numel(CSidx.(fn{ii}))
        x = GTS.UpperBoundary(CSidx.(fn{ii})(jj));
        w = GTS.LowerBoundary(CSidx.(fn{ii})(jj))-GTS.UpperBoundary(CSidx.(fn{ii})(jj));
        rectangle('Position',[x,y,w,h],'EdgeColor','none','FaceColor',cm(ii,:))
    end
end
ylim([8,42])
geologictimescale(0,GTS.LowerBoundary(DArange(end)),...
    'normal','reverse',gca,'standard','stages','off',4,1)  
ax1 = gca; ax1.FontSize = 11; ax1.FontName = 'helvetica'; ax1.YTick = [];
ylabel('Climate state','FontName','helvetica','FontSize',13,'FontWeight','bold')
xlabel('Age (Ma)','FontName','helvetica','FontSize',13,'FontWeight','bold')
text(482,38,'A','FontName','helvetica','FontSize',15,'FontWeight','bold','color','w')

% (c) Distributions
ax2 = subplot(3,2,[3,5]);
ax2.Position = [.06,.075,.425,.4];
hold on
xedge = [floor(min(mgmst))-8:.25:ceil(max(mgmst))+8]';
y = ksdensity(mgmst,xedge);
xtext = [12,16,24,32.75,37];
ytext = [.0375,.05325,.069,.05325,.0375]';
statelab = ["Coldhouse","Coolhouse","Transitional","Warmhouse","Hothouse"];
for ii = 1:numel(fn)
    if ii == 1
        idx = xedge <= P(ii);
        temptext = sprintf("(%.0f-%.0f%sC)",round(min(mgmst)),round(P(ii)),char(176));
    elseif ii == numel(fn)
        idx = xedge > P(ii-1);
        temptext = sprintf("(%.0f-%.0f%sC)",round(P(ii-1)),round(max(mgmst)),char(176));
    else
        idx = xedge > P(ii-1) & xedge <= P(ii);
        temptext = sprintf("(%.0f-%.0f%sC)",round(P(ii-1)),round(P(ii)),char(176));
    end
    plotx = [xedge(idx);flipud(xedge(idx))];
    ploty = [y(idx);zeros(sum(idx),1)];
    fill(plotx,ploty,cm(ii,:),'FaceAlpha',.75,'EdgeColor','none');
    text(xtext(ii),ytext(ii),statelab(ii),'FontWeight','bold','FontName',...
        'helvetica','FontSize',11,'HorizontalAlignment','center','Color',cm_font(ii,:))
    text(xtext(ii),ytext(ii)-.0045,temptext,'FontName','helvetica','FontSize',...
        11,'HorizontalAlignment','center','Color',cm_font(ii,:))
end
ylim([0 .069]),xlim([3 44])
ax2.FontName = 'helvetica'; ax2.FontSize = 11;
xlabel(['GMST (',char(176),'C)'],'FontName','helvetica','FontSize',13,...
    'FontWeight','bold')
ylabel('Density','FontName','helvetica','FontSize',13,'FontWeight','bold')
text(3.75,.0665,'B','FontName','helvetica','FontSize',15,'FontWeight','bold')
% (d) Pie chart
ax3b = axes('Position',[.6,.2,.1,.1])
ax3b.Position = [.3,.425,.2,.2];
h = pie(CStime.Total(1,:));
for ii = 2:2:2*numel(CStime.Total(1,:))
    h(ii).Position = 0.4*h(ii).Position;
    if ii ~= 6
        h(ii).Color = 'w';
    end
    h(ii).FontName = 'helvetica';
    h(ii).FontWeight = 'bold';
    h(ii).FontSize = 11;
    h(ii-1).EdgeColor = 'w';
end
colormap(cm)
% manually adjust labels to be in center of each slice
h(2).Position = [-.005,.475,0];
h(4).Position = [-.35,.1,0];
h(8).Position = [.3,-.175,0];
h(10).Position = [.125,.4,0];

% (d) Plot latitudinal gradients
ax4 = subplot(3,6,[10:12,16:18]);
ax4.Position = [.575,.075,.4,.55];
hold on, box on
for ii = 1:numel(fn)
    % 1 sigma bound
    y=[];
    for jj = 1:numel(CSidx.(fn{ii}))
        y = [y, LTG{CSidx.(fn{ii})(jj)}];
    end
    f = fill([Lat;flipud(Lat)],[prctile(y,16,2);...
        flipud(prctile(y,84,2))],cm(ii,:),...
        'FaceAlpha',.25,'EdgeColor','none');
    plot(Lat,median(y,2),'-','color',cm_font(ii,:),'LineWidth',3)
end
% Tidy plot
xlim([-90,90]);
ax4.XTick = [-90:30:90];
ax4.FontName = 'helvetica';
ax4.FontSize = 11;
xlabel('Latitude','FontName','helvetica','FontSize',13,'FontWeight','bold')
ylabel(['Air temperautre (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')
ylim([-50 46])
text(-86,42,'C','FontName','helvetica','FontSize',15,'FontWeight','bold','color','k')
ax4.Color = 'none';
% Move pie chart & Tidy
ax3b.Position = [.675,.1,.2,.25];
text(ax3b, 0,-1, sprintf('Proportion of time\nspent in each state'),...
    'FontName','helvetica','FontSize',11,'FontWeight','normal', ...
    'HorizontalAlignment','center','VerticalAlignment','top');
box(ax3b,'on')
ax3b.Visible = 'on';ax3b.XTick = [];ax3b.YTick = [];
ylim(ax3b,[-1.7,1.1])
xlim(ax3b,[-1.3,1.3])
text(ax3b,-1.2,.8,'D','FontName','helvetica','FontSize',15,'FontWeight','bold')

% (e) Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end
