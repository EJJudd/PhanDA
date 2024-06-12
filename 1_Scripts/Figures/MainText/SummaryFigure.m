

% PART 1: LOAD DATA
% Directory details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
assdate = '27Jul2023';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","Ndata","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")
load("GTS2020_PETM.mat")

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
curveprctiles = [5,10,16,25,40,50,60,75,84,90,95];
pgmst = cell2mat(cellfun(@(x) prctile(x, curveprctiles), GMST, 'UniformOutput', false));
% subdivide GMST by climate state
stateprctiles = [0:20:100];
P = prctile(pgmst(:,curveprctiles==50), stateprctiles);
CS = discretize(pgmst(:,curveprctiles==50),P);


%% Part 3: MAKE FIGURE
cm = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);
fig=figure('Name','SummaryGMST','color','w','Units','centimeters'); 
fig.Position=[-8 14 33 11];
pause(0.5)
% (b) Plot GMST
ax = gca; hold on
% Smooth curve (for broader audiance)
x = []; y = []; 
y5 = []; y10 = []; y16 = []; y25 = []; y40 = [];
y60 = []; y75 = []; y84 = []; y90 = []; y95 = [];
for ii = 1:numel(GTS.Average)-1
    if ii == 46 || ii == 63 || ii == 21 || ii == 49
    x = [x;linspace(GTS.Average(ii),GTS.Average(ii+1),5)'];
    y = [y;linspace(pgmst(ii,curveprctiles==50),pgmst(ii+1,curveprctiles==50),5)'];
    y5 = [y5;linspace(pgmst(ii,curveprctiles==5),pgmst(ii+1,curveprctiles==5),5)'];
    y10 = [y10;linspace(pgmst(ii,curveprctiles==10),pgmst(ii+1,curveprctiles==10),5)'];
    y16 = [y16;linspace(pgmst(ii,curveprctiles==16),pgmst(ii+1,curveprctiles==16),5)'];
    y25 = [y25;linspace(pgmst(ii,curveprctiles==25),pgmst(ii+1,curveprctiles==25),5)'];
    y40 = [y40;linspace(pgmst(ii,curveprctiles==40),pgmst(ii+1,curveprctiles==40),5)'];
    y60 = [y60;linspace(pgmst(ii,curveprctiles==60),pgmst(ii+1,curveprctiles==60),5)'];
    y75 = [y75;linspace(pgmst(ii,curveprctiles==75),pgmst(ii+1,curveprctiles==75),5)'];
    y84 = [y84;linspace(pgmst(ii,curveprctiles==84),pgmst(ii+1,curveprctiles==84),5)'];
    y90 = [y90;linspace(pgmst(ii,curveprctiles==90),pgmst(ii+1,curveprctiles==90),5)'];
    y95 = [y95;linspace(pgmst(ii,curveprctiles==95),pgmst(ii+1,curveprctiles==95),5)'];
    else
    x = [x;linspace(GTS.Average(ii),GTS.Average(ii+1),3)'];
    y = [y;linspace(pgmst(ii,curveprctiles==50),pgmst(ii+1,curveprctiles==50),3)'];
    y5 = [y5;linspace(pgmst(ii,curveprctiles==5),pgmst(ii+1,curveprctiles==5),3)'];
    y10 = [y10;linspace(pgmst(ii,curveprctiles==10),pgmst(ii+1,curveprctiles==10),3)'];
    y16 = [y16;linspace(pgmst(ii,curveprctiles==16),pgmst(ii+1,curveprctiles==16),3)'];
    y25 = [y25;linspace(pgmst(ii,curveprctiles==25),pgmst(ii+1,curveprctiles==25),3)'];
    y40 = [y40;linspace(pgmst(ii,curveprctiles==40),pgmst(ii+1,curveprctiles==40),3)'];
    y60 = [y60;linspace(pgmst(ii,curveprctiles==60),pgmst(ii+1,curveprctiles==60),3)'];
    y75 = [y75;linspace(pgmst(ii,curveprctiles==75),pgmst(ii+1,curveprctiles==75),3)'];
    y84 = [y84;linspace(pgmst(ii,curveprctiles==84),pgmst(ii+1,curveprctiles==84),3)'];
    y90 = [y90;linspace(pgmst(ii,curveprctiles==90),pgmst(ii+1,curveprctiles==90),3)'];
    y95 = [y95;linspace(pgmst(ii,curveprctiles==95),pgmst(ii+1,curveprctiles==95),3)'];
    end
end
x = unique(x,'stable');
y = unique(y,'stable');
y5 = unique(y5,'stable');
y10 = unique(y10,'stable');
y16 = unique(y16,'stable');
y25 = unique(y25,'stable');
y40 = unique(y40,'stable');
y60 = unique(y60,'stable');
y75 = unique(y75,'stable');
y84 = unique(y84,'stable');
y90 = unique(y90,'stable');
y95 = unique(y95,'stable');
xx = linspace(0, round(GTS.Average(end)), round(GTS.Average(end) * 100))';
yy = spline(x,y,xx);
yy5 = spline(x,y5,xx);
yy10 = spline(x,y10,xx);
yy16 = spline(x,y16,xx);
yy25 = spline(x,y25,xx);
yy40 = spline(x,y40,xx);
yy60 = spline(x,y60,xx);
yy75 = spline(x,y75,xx);
yy84 = spline(x,y84,xx);
yy90 = spline(x,y90,xx);
yy95 = spline(x,y95,xx);
xx(end+1) = GTS.LowerBoundary(end);
yy(end+1) = yy(end);
yy5(end+1) = yy5(end);
yy10(end+1) = yy10(end);
yy16(end+1) = yy16(end);
yy25(end+1) = yy25(end);
yy40(end+1) = yy40(end);
yy60(end+1) = yy60(end);
yy75(end+1) = yy75(end);
yy84(end+1) = yy84(end);
yy90(end+1) = yy90(end);
yy95(end+1) = yy95(end);

fill([xx;flipud(xx)],[yy5;flipud(yy95)],...
    'w','FaceColor',[.9,.9,.9],'EdgeColor','none');
fill([xx;flipud(xx)],[yy10;flipud(yy90)],...
    'w','FaceColor',[.8,.8,.8],'EdgeColor','none');
fill([xx;flipud(xx)],[yy16;flipud(yy84)],...
    'w','FaceColor',[.7,.7,.7],'EdgeColor','none');
fill([xx;flipud(xx)],[yy25;flipud(yy75)],...
    'w','FaceColor',[.6,.6,.6],'EdgeColor','none');
fill([xx;flipud(xx)],[yy40;flipud(yy60)],...
    'w','FaceColor',[.5,.5,.5],'EdgeColor','none');
plot(xx,yy,'k-','LineWidth',2.5)
ax.YTick = [10:5:40];

% Add climate states
y = 44.5; h = 3.5;
for ii = 1:numel(CS)
    x = GTS.UpperBoundary(ii);
    w = GTS.LowerBoundary(ii) - x;
    rectangle('Position',[x,y,w,h],'FaceColor',cm(CS(ii),:),'EdgeColor','none')
end
rectangle('Position',[0,y,GTS.LowerBoundary(end),h],'EdgeColor',[0,0,0],'FaceColor','none','LineWidth',1.5)

% Tidy plot and add time scale
set(gca,'XDir','reverse','FontName','Arial','FontSize',11)
xlim([0,GTS.LowerBoundary(end)])
ylim([3 48])
rectangle('Position',[0,3,GTS.LowerBoundary(22),3.5],'FaceColor',[.25,.25,.25],'EdgeColor','k','LineWidth',1.5)
rectangle('Position',[GTS.LowerBoundary(22),3,GTS.LowerBoundary(50)-GTS.LowerBoundary(22),3.5],'FaceColor',[.35,.35,.35],'EdgeColor','k','LineWidth',1.5)
rectangle('Position',[GTS.LowerBoundary(50),3,GTS.LowerBoundary(end)-GTS.LowerBoundary(50),3.5],'FaceColor',[.45,.45,.45],'EdgeColor','k','LineWidth',1.5)
text(mean([0,GTS.LowerBoundary(22)]),4.75,'Cenozoic','FontName','Arial',...
    'FontSize',13,'Color','w','HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontWeight','bold')
text(mean([GTS.LowerBoundary(22),GTS.LowerBoundary(50)]),4.75,'Mesozoic','FontName','Arial',...
    'FontSize',13,'Color','w','HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontWeight','bold')
text(mean([GTS.LowerBoundary(50),GTS.LowerBoundary(end)]),4.75,'Paleozoic','FontName','Arial',...
    'FontSize',13,'Color','w','HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontWeight','bold')
ax.XAxis.MinorTick = 'on';
ax.XAxis.TickDirection = 'out';
xlabel('Age (millions of years ago)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel(['Global temperature (',char(176),'C)'],'FontName','Arial','FontSize',15,'FontWeight','bold')

% Add legends
x = 280; w = 60;
y1 = 35; y2 = 38; h = 2;
for ii = 1:5
    rectangle('Position',[x+ii*(w/5)-(w/5),y1,w/5,h],...
        'FaceColor',[.4,.4,.4]+ii/10,'EdgeColor','none')
    rectangle('Position',[x+ii*(w/5)-(w/5),y2,w/5,h],...
        'FaceColor',cm(end-ii+1,:),'EdgeColor','none')
end
rectangle('Position',[x,y1,w,h],'EdgeColor',[0,0,0],'FaceColor','none','LineWidth',1)
text(x+w,y1-.5,'less likely ','FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold')
text(x,y1-.5,'more likely ','FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold')

rectangle('Position',[x,y2,w,h],'EdgeColor',[0,0,0],'FaceColor','none','LineWidth',1)
text(x+w,y2+h+.5,'cooler climate ','FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold')
text(x,y2+h+.5,'warmer climate ','FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold')


rectangle('Position',[0,3,GTS.LowerBoundary(end),45],'EdgeColor',[0,0,0],'FaceColor','none','LineWidth',1.5)
annotation('textbox',[.125,.945,0,0],'String',['Climate',newline,'state'],'FontName','Arial','FontSize',15,'FontWeight','bold','HorizontalAlignment','right')

%export_fig(gcf,[figdir,'/SummaryFigure.tiff'],'-p0.01','-m5')


%% GRAVEYARD

% Add color stripes
deltagmst = cell2mat(cellfun(@(x) median(x)-median(GMST{1}), GMST, 'UniformOutput', false));
deltabins = [linspace(-2.75,0,9),linspace(0,22,9)];
deltabins(10) = [];
deltadisc = discretize(deltagmst,deltabins);

cm = hex2rgb({'#08306b', '#08519c', '#2171b5', '#4292c6','#6baed6', ...
    '#9ecae1', '#c6dbef', '#deebf7','#fee0d2', '#fcbba1', '#fc9272', ...
    '#fb6a4a','#ef3b2c', '#cb181d', '#a50f15', '#67000d'},1);
y = 5; h = 40;
for ii = 1:numel(deltagmst)
    x = GTS.UpperBoundary(ii);
    w = GTS.LowerBoundary(ii) - x;
    rectangle('Position',[x,y,w,h],'FaceColor',cm(deltadisc(ii),:),'EdgeColor','none')
end
% Add colormap
%rectangle('Position',[223,36.4,144,43.9],'FaceColor',[.9,.9,.9,.6],'EdgeColor','none')
rectangle('Position',[13,7.71,144,7.59],'FaceColor',[.9,.9,.9,.6],'EdgeColor','none')
cmbins = round(deltabins,1);
cmplot = zeros(range(cmbins)*10,3);
xcm = [cmbins(1):.1:cmbins(end)];
xcmbin = discretize(xcm,cmbins);
cmplot = cm(xcmbin,:);
c = colorbar('south');
pause(1)
colormap(cmplot)
c.Position = [.667,.22,.21,.04];
%c.Position = [.33,.87,.21,.04];
set(c,'Color','w','FontName','Arial','FontSize',11,'TickDirection','out','LineWidth',1.5)
caxis([-2.8,22])
c.Ticks = [cmbins(1),cmbins([9:2:end])];
ylabel(c,['\DeltaGMST (',char(176),'C relative to preindustrial)'],'FontName','Arial','FontSize',13,'FontWeight','bold','Color','w')
