%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GMST vs CO2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          
% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
figname = 'Fig4_CO2&Forcings.png';
savefig = true;        

% Select colormap
cm_grey = customcolormap(linspace(0,1,3),{'#515151','#B8B8B8','#ECECEC'},5);
cm_col = hex2rgb(['#0a9396';'#ffb703';'#ca6702'],1);

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Documents/PhanDA/5_Outputs/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","Ndata","Index","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")
load("GTS2020_PETM.mat")
load("PhanerozoicCO2v9.mat", "PhanerozoicCO2")

% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
swCorr = ["snowball","off"];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",swCorr)) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);
% Revise CO2 & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
PhanerozoicCO2 = combinestages(PhanerozoicCO2,"CO2",Preferences,endsize,2);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
% Calculate percentiles
P = [5,10,16,25,40,50,60,75,84,90,95];
Pgmst = cell2mat(cellfun(@(x) prctile(x, P), GMST, 'UniformOutput', false));
CO2 = prctile(PhanerozoicCO2,P,2);

% Define eras
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
mesozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Selandian/Danian") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Induan"));
paleozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Tremadocian"));
paleonocopse = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));
allnocopse = find(GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));

%% PART 3: MAKE FIGURE
% (a) Initialize figure
fig = figure('Units','inches','Position',[1,1,6.5*1.5,5.75*1.5],'Color','w');
tiledlayout(6,2,'Padding','none','TileSpacing','compact');
fig.Units = 'pixels';
% (b) Plot GMST & CO2
ax1 = nexttile([3,2]); hold on
for ii = 1:5
    fill([GTS.Average;flipud(GTS.Average)],...
        [Pgmst(:,ii);flipud(Pgmst(:,end+1-ii))],'k','FaceColor',cm_grey(ii,:),'EdgeColor','none');
end
plot(GTS.Average,cell2mat(cellfun(@(x) median(x), GMST, 'UniformOutput', false)),'k-','LineWidth',2)
ylim([-42,44.9])
ax1.YTick = [10:10:40];
ylabel(['                 GMST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold','HorizontalAlignment','left');
% (c) Plot CO2
yyaxis right
for ii = 1:5
    fill([GTS.Average(allnocopse);flipud(GTS.Average(allnocopse))],...
        log([CO2(allnocopse,ii);flipud(CO2(allnocopse,end+1-ii))]),'k','FaceColor',cm_grey(ii,:),'EdgeColor','none');
end
plot(GTS.Average(allnocopse),log(median(PhanerozoicCO2(allnocopse,:),2)),'k-','LineWidth',2)
ylim(log([25,100000]))
ax1.YTick = log([100,200,500,1000,2000,5000]);
ax1.YAxis(2).MinorTickValues = log([100:100:1000,2000:1000:5000]);
ax1.YTickLabel = ["100","200","500","1000","2000","5000"];
ax1.FontSize = 11; 
ax1.FontName = 'helvetica';
ylabel(['CO','\fontsize{7}2', '\fontsize{13}  (ppmv)'],'FontName','helvetica','FontWeight','bold','Interpreter','tex')
ax1.YColor = 'k';
caxis([5 95])
yyaxis left
geologictimescale(0,GTS.LowerBoundary(end),...
    'normal','reverse',gca,'standard','stages','off',7.5,2,'helvetica',11,'k','k')
xlabel('Age (Ma)','FontName','helvetica','FontSize',13,'FontWeight','bold')
% text(ax1,15,40,'A','FontName','helvetica','FontSize',15,'FontWeight','bold','color','k')
% (d) Plot GMST vs. CO2
ax2 = nexttile([3,1]); hold on, box on
errorbar(log2(CO2(allnocopse,P(1,:)==50)),...
        Pgmst(allnocopse,P(1,:)==50),...
        Pgmst(allnocopse,P(1,:)==50) - Pgmst(allnocopse,P(1,:)==16),...
        Pgmst(allnocopse,P(1,:)==84) - Pgmst(allnocopse,P(1,:)==50),...
        log2(CO2(allnocopse,P(1,:)==84)) - log2(CO2(allnocopse,P(1,:)==50)),...
        log2(CO2(allnocopse,P(1,:)==50)) - log2(CO2(allnocopse,P(1,:)==16)), ...
        '.','Color',[.65,.65,.65],'LineWidth',1.5,'CapSize',0)
p2 = plot(log2(CO2(mesozoic,P(1,:)==50)),...
    Pgmst(mesozoic,P(1,:)==50),'ko',...
    'MarkerFaceColor',hex2rgb('#FFCF57',1),'MarkerSize',10);
p1 = plot(log2(CO2(cenozoic,P(1,:)==50)),Pgmst(cenozoic,P(1,:)==50),'ko',...
    'MarkerFaceColor',hex2rgb('#5CB7B9',1),'MarkerSize',10);
p3 = plot(log2(CO2(paleonocopse,P(1,:)==50)),Pgmst(paleonocopse,P(1,:)==50),'ko',...
    'MarkerFaceColor',hex2rgb('#DC9A56',1),'MarkerSize',10);
% York regression
x = log2(PhanerozoicCO2);
y = cell2mat(cellfun(@(x) datasample(x, numel(GMST{1}), 'Replace', false)', ...
    GMST, 'UniformOutput', false));
X = [min(median(x(allnocopse,:),2)),max(median(x(allnocopse,:),2))];
r = 0;
% All
sx = range(prctile(x(allnocopse,:),[16,84],2),2)./2;
sy = range(prctile(y(allnocopse,:),[16,84],2),2)./2;
[b, m, ~, sm] = york_fit(median(x(allnocopse,:),2)',median(y(allnocopse,:),2)',sx',sy',r);
plot(X,X*m+b,'k--','LineWidth',2)
text(log2(175),38,'Phanerozoic','FontName','helvetica','FontSize',13,'FontWeight','bold','color','k')
text(log2(175),38-1.75,['AESS = ',sprintf('%0.1f',m),'\pm',sprintf('%0.1f%sC',sm,char(176))],...
    'FontName','helvetica','FontSize',11,'FontWeight','bold','color','k')
text(log2(175),38-1.75*2,sprintf('r = %.02f',corr(median(x(allnocopse,:),2),...
    median(y(allnocopse,:),2))),'FontName','helvetica','FontSize',11,...
    'FontWeight','bold','color','k')
% Annotate
text(log2(1500),15.25,"Cenozoic",'FontName','helvetica','FontSize',13,'Color',cm_col(1,:),'FontWeight','bold')
text(log2(1500),13.5,"Mesozoic",'FontName','helvetica','FontSize',13,'Color',cm_col(2,:),'FontWeight','bold')
text(log2(1500),11.75,"Paleozoic",'FontName','helvetica','FontSize',13,'Color',cm_col(3,:),'FontWeight','bold')
text(log2(175),41,'B','FontName','helvetica','FontSize',15,'FontWeight','bold','color','k')
% Tidy
ylim([8 43])
xlim(ax2,log2([160 4500]))
ax2.XTick = log2([100,200,500,1000,2000,5000]);
ax2.XTickLabel = [100,200,500,1000,2000,5000];
ax2.XMinorTick = 'on';
ax2.XAxis.MinorTickValues = log2([100:100:1900,2000:1000:5000]);
ax2.TickLength = [.0225,.05];
ax2.FontSize = 11; 
ax2.FontName = 'helvetica';
xlabel(['\fontsize{13}CO','\fontsize{7}2', '\fontsize{13}  (ppmv)'],'FontName','helvetica','FontWeight','bold','Interpreter','tex')
ylabel(['GMST (',char(176),'C)'],'FontName','helvetica','FontSize',13,'FontWeight','bold')

% (e) Plot CO2 by climate state
% subdivide GMST by percentile
ax3 = nexttile([3,1]); hold on, box on
P2 = [20:20:100];
mgmst = cell2mat(cellfun(@(x) nanmedian(x), GMST, 'UniformOutput', false));
P2 = prctile(mgmst, P2);
cm_grey = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);
cm_font = cm_grey; cm_font(3,:) = [.5,.5,.5];
CSidx.ih = find(mgmst(allnocopse)<=P2(1));
CSidx.ch = find(mgmst(allnocopse)>P2(1) & mgmst(allnocopse)<=P2(2));
CSidx.tr = find(mgmst(allnocopse)>P2(2) & mgmst(allnocopse)<=P2(3));
CSidx.gh = find(mgmst(allnocopse)>P2(3) & mgmst(allnocopse)<=P2(4));
CSidx.hh = find(mgmst(allnocopse)>P2(4));
cm_grey = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);
fn = fieldnames(CSidx);
w = 1;
lab = ["coldhouse","coolhouse","transitional","warmhouse","hothouse"];
sf = linspace(.3,.9,5);
allnocopse = find(GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));
for ii = 1:numel(fn)
    x = ii;
    y = prctile(PhanerozoicCO2(intersect(CSidx.(fn{ii}),allnocopse),:),[5 16 50 84 95],[2,1]);   
    ym = prctile(PhanerozoicCO2(setdiff(intersect(CSidx.(fn{ii}),allnocopse),mesozoic),:),[5 16 50 84 95],[2,1]);   
    rectangle('Position',[x,y(1),w,y(end)-y(1)],'FaceColor',[cm_grey(ii,:),.25],'EdgeColor','none')
    rectangle('Position',[x,y(2),w,y(end-1)-y(1)],'FaceColor',[cm_grey(ii,:),.35],'EdgeColor','none')
    plot([x,x+w],[y(3),y(3)],'-','Color',cm_grey(ii,:),'LineWidth',3)
    plot([x,x+w],[ym(3),ym(3)],'--','Color',cm_grey(ii,:),'LineWidth',1.5)
    if ii<4
        text(ii+.5*w,y(3)+100*sf(ii),sprintf('%.0f',y(3)),'FontName','helvetica','FontSize',10,'Color',cm_font(ii,:),'HorizontalAlignment','center')
        text(ii+.5*w,ym(3)-100*sf(ii),sprintf('%.0f',ym(3)),'FontName','helvetica','FontSize',10,'Color',cm_font(ii,:),'HorizontalAlignment','center')
    else
        text(ii+.5*w,y(3)-100*sf(ii),sprintf('%.0f',y(3)),'FontName','helvetica','FontSize',10,'Color',cm_font(ii,:),'HorizontalAlignment','center')
        text(ii+.5*w,ym(3)+100*sf(ii),sprintf('%.0f',ym(3)),'FontName','helvetica','FontSize',10,'Color',cm_font(ii,:),'HorizontalAlignment','center')
    end        
    text(ii+.5*w,y(1)-75*sf(ii),lab(ii),'FontWeight','bold','FontName','helvetica','FontSize',10,'Color',cm_font(ii,:),'HorizontalAlignment','center')
    fprintf('%s = %.0f ppm (%.0f - %.0f)\n',fn{ii},y(3),y(2),y(4))
end
set(ax3,'YScale','log')
ylim([150,3250])
ax3.YTick = [200,500,1000,2000];
ax3.YMinorTick = 'on';
ax3.YAxis.MinorTickValues = [200:100:4000];    
ax3.TickLength = [.0225,.05];
xlim([.9 6.1])
ax3.XTick = [];
ax3.FontSize = 11;
ax3.FontName = 'helvetica';
ylabel('CO_2 (ppmv)','FontName','helvetica','FontSize',13,'FontWeight','bold')
xlabel('Climate state','FontName','helvetica','FontSize',13,'FontWeight','bold')
text(ax3,1.05,2800,'C','FontName','helvetica','FontSize',15,'FontWeight','bold','color','k')

if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end

