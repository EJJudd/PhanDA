%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SupFig  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          
% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
figname = 'SupFig_ACS.png';
savefig = true;        

% Select colormap
cm = hex2rgb({'#004F60';'#0a9396';'#ffb703';'#ca6702';'#9b2226'},1);

% PART 1: LOAD DATA
% Directory details
assdate = '27Jul2023';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
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
sbCorr = [true, false];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",string(pHCorr))) & ...
    contains(ItName,strcat("SnowballCorr = ",string(sbCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);
% Revise CO2 & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
PhanerozoicCO2 = combinestages(PhanerozoicCO2,"CO2",Preferences,endsize,2);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
% Calculate percentiles
P = [5:1:95];
Pgmst = cell2mat(cellfun(@(x) prctile(x, P), GMST, 'UniformOutput', false));
CO2 = prctile(PhanerozoicCO2,P,2);
[P,age_grid] = meshgrid(P,GTS.Average);

% Define eras
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
mesozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Selandian/Danian") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Induan"));
paleozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Tremadocian"));
paleonocopse = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));
allnocopse = find(GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));
allnocopsemeso = setdiff(find(GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian")),mesozoic);

%% PART 3: MAKE FIGURE
% (a) Initialize figure
fig = figure('Units','inches','Position',[1,1,12.5,4.125],'Color','w');
tiledlayout(1,3,'Padding','none','TileSpacing','compact');

ax1 = nexttile; hold on, box on
errorbar(log2(CO2(cenozoic,P(1,:)==50)),...
        Pgmst(cenozoic,P(1,:)==50),...
        Pgmst(cenozoic,P(1,:)==50) - Pgmst(cenozoic,P(1,:)==16),...
        Pgmst(cenozoic,P(1,:)==84) - Pgmst(cenozoic,P(1,:)==50),...
        log2(CO2(cenozoic,P(1,:)==84)) - log2(CO2(cenozoic,P(1,:)==50)),...
        log2(CO2(cenozoic,P(1,:)==50)) - log2(CO2(cenozoic,P(1,:)==16)), ...
        '.','Color',[.65,.65,.65],'LineWidth',2,'CapSize',0)
plot(log2(CO2(cenozoic,P(1,:)==50)),Pgmst(cenozoic,P(1,:)==50),'ko',...
    'MarkerFaceColor',hex2rgb('#66BBBD',1),'MarkerSize',10);
% York regression
X = log2(CO2(:,P(1,:) == 50));
Y = Pgmst(:,P(1,:) == 50);
sx = (log2(CO2(cenozoic,P(1,:)==84)) - log2(CO2(cenozoic,P(1,:)==16)))./2;
sy = (Pgmst(cenozoic,P(1,:)==84) - Pgmst(cenozoic,P(1,:)==16))./2;
r = 0;
[b, m, ~, sm] = york_fit(X(cenozoic)',...
    Y(cenozoic)',sx',sy', r);
x = min(X(cenozoic)-sx):.1:max(X(cenozoic)+sx);
y = m.*x + b;
plot(x,y,'k--','LineWidth',2)
ylim([8 41])
xlim(log2([180 4500]))
text(log2(195),39.25,'A','FontName','Arial','FontSize',15,'FontWeight','bold','color','k')
text(log2(195),37,'Cenozoic','FontName','Arial','FontSize',13,'FontWeight','bold','color',cm(2,:))
text(log2(195),35.5,['ACS = ',sprintf('%0.1f',m),'\pm',sprintf('%0.1f%sC',sm,char(176))],...
    'FontName','Arial','FontSize',11,'FontWeight','bold','color',cm(2,:))
text(log2(195),34.25,['      r = ',sprintf('%.02f',corr(X(cenozoic),Y(cenozoic)))],...
    'FontName','Arial','FontSize',11,'FontWeight','bold','color',cm(2,:))
ax1.XTick = log2([100,200,500,1000,2000,5000]);
ax1.XTickLabel = [100,200,500,1000,2000,5000];
ax1.XMinorTick = 'on';
ax1.XAxis.MinorTickValues = log2([100:100:1900,2000:1000:5000]);
ax1.TickLength = [.0225,.05];
xlabel(['\fontsize{13}CO','\fontsize{7}2', '\fontsize{13}  (ppmv)'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')

% Paleozoic
ax1 = nexttile; hold on, box on
errorbar(log2(CO2(paleonocopse,P(1,:)==50)),...
        Pgmst(paleonocopse,P(1,:)==50),...
        Pgmst(paleonocopse,P(1,:)==50) - Pgmst(paleonocopse,P(1,:)==16),...
        Pgmst(paleonocopse,P(1,:)==84) - Pgmst(paleonocopse,P(1,:)==50),...
        log2(CO2(paleonocopse,P(1,:)==84)) - log2(CO2(paleonocopse,P(1,:)==50)),...
        log2(CO2(paleonocopse,P(1,:)==50)) - log2(CO2(paleonocopse,P(1,:)==16)), ...
        '.','Color',[.65,.65,.65],'LineWidth',2,'CapSize',0)
plot(log2(CO2(paleonocopse,P(1,:)==50)),Pgmst(paleonocopse,P(1,:)==50),'ko',...
    'MarkerFaceColor',hex2rgb('#DEA061',1),'MarkerSize',10);
% York regression
sx = (log2(CO2(paleonocopse,P(1,:)==84)) - log2(CO2(paleonocopse,P(1,:)==16)))./2;
sy = (Pgmst(paleonocopse,P(1,:)==84) - Pgmst(paleonocopse,P(1,:)==16))./2;
r = 0;
[b, m, ~, sm] = york_fit(X(paleonocopse)',...
    Y(paleonocopse)',sx',sy', r);
x = min(X(paleonocopse)-sx/2):.1:max(X(paleonocopse)+sx/2);
y = m.*x + b;
plot(x,y,'k--','LineWidth',2)
fprintf('Paleozoic ACS = %.1f %s %.1f \n', m, char(177), sm)
ylim([8 41])
xlim(log2([180 4500]))
text(log2(195),39.25,'B','FontName','Arial','FontSize',15,'FontWeight','bold','color','k')
text(log2(195),37,'Paleozoic','FontName','Arial','FontSize',13,'FontWeight','bold','color',cm(4,:))
text(log2(195),35.5,['ACS = ',sprintf('%0.1f',m),'\pm',sprintf('%0.1f%sC',sm,char(176))],...
    'FontName','Arial','FontSize',11,'FontWeight','bold','color',cm(4,:))
text(log2(195),34.25,['      r = ',sprintf('%.02f',corr(X(paleonocopse),Y(paleonocopse)))],...
    'FontName','Arial','FontSize',11,'FontWeight','bold','color',cm(4,:))
ax1.XTick = log2([100,200,500,1000,2000,5000]);
ax1.XTickLabel = [100,200,500,1000,2000,5000];
ax1.XMinorTick = 'on';
ax1.XAxis.MinorTickValues = log2([100:100:1900,2000:1000:5000]);
ax1.TickLength = [.0225,.05];
xlabel(['\fontsize{13}CO','\fontsize{7}2', '\fontsize{13}  (ppmv)'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')

% Paleozoic & Cenozoic
ax3 = nexttile; hold on, box on
errorbar(log2(CO2(allnocopsemeso,P(1,:)==50)),...
        Pgmst(allnocopsemeso,P(1,:)==50),...
        Pgmst(allnocopsemeso,P(1,:)==50) - Pgmst(allnocopsemeso,P(1,:)==16),...
        Pgmst(allnocopsemeso,P(1,:)==84) - Pgmst(allnocopsemeso,P(1,:)==50),...
        log2(CO2(allnocopsemeso,P(1,:)==84)) - log2(CO2(allnocopsemeso,P(1,:)==50)),...
        log2(CO2(allnocopsemeso,P(1,:)==50)) - log2(CO2(allnocopsemeso,P(1,:)==16)), ...
        '.','Color',[.65,.65,.65],'LineWidth',2,'CapSize',0)
plot(log2(CO2(cenozoic,P(1,:)==50)),Pgmst(cenozoic,P(1,:)==50),'ko',...
    'MarkerFaceColor',hex2rgb('#66BBBD',1),'MarkerSize',10);
plot(log2(CO2(paleonocopse,P(1,:)==50)),Pgmst(paleonocopse,P(1,:)==50),'ko',...
    'MarkerFaceColor',hex2rgb('#DEA061',1),'MarkerSize',10);
% York regression
sx = (log2(CO2(allnocopsemeso,P(1,:)==84)) - log2(CO2(allnocopsemeso,P(1,:)==16)))./2;
sy = (Pgmst(allnocopsemeso,P(1,:)==84) - Pgmst(allnocopsemeso,P(1,:)==16))./2;
r = 0;
[b, m, ~, sm] = york_fit(X(allnocopsemeso)',...
    Y(allnocopsemeso)',sx',sy', r);
x = min(X(allnocopsemeso)-sx/2):.1:max(X(allnocopsemeso)+sx/2);
y = m.*x + b;
plot(x,y,'k--','LineWidth',2)
ylim([8 41])
xlim(log2([180 4500]))
text(log2(195),39.25,'C','FontName','Arial','FontSize',15,'FontWeight','bold','color','k')
text(log2(195),37,'Cenozoic & Paleozoic','FontName','Arial','FontSize',13,'FontWeight','bold','color','k')
text(log2(195),35.5,['ACS = ',sprintf('%0.1f',m),'\pm',sprintf('%0.1f%sC',sm,char(176))],...
    'FontName','Arial','FontSize',11,'FontWeight','bold','color','k')
text(log2(195),34.25,['      r = ',sprintf('%.02f',corr(X(allnocopsemeso),Y(allnocopsemeso)))],...
    'FontName','Arial','FontSize',11,'FontWeight','bold','color','k')
ax3.XTick = log2([100,200,500,1000,2000,5000]);
ax3.XTickLabel = [100,200,500,1000,2000,5000];
ax3.XMinorTick = 'on';
ax3.XAxis.MinorTickValues = log2([100:100:1900,2000:1000:5000]);
ax3.TickLength = [.0225,.05];
xlabel(['\fontsize{13}CO','\fontsize{7}2', '\fontsize{13}  (ppmv)'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')

if savefig
    export_fig(gcf,[figdir,'/Supplemental/',figname],'-p0.01','-m5')
end

