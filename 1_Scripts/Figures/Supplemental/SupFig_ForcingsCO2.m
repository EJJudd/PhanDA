%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Forcings    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          
% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
savefig = false;        

% Select colormap
cm_grey = customcolormap(linspace(0,1,2),{'#6A6A6A','#ECECEC'},50);
cm_grey = [cm_grey;cm_grey(end,:);flipud(cm_grey)];
cm_col = hex2rgb(['#0a9396';'#ffb703';'#ca6702'],1);

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Documents/PhanDA/5_Outputs/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","Ndata","Index","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("PhanerozoicCO2v9.mat", "PhanerozoicCO2")
load("ExpNo.mat","ExpNo")
load("lsmask.mat","lsmask");
lsmask = struct2cell(lsmask);


% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
swCorr = ["snowball", "veizer","off"];
%swCorr = "veizer";
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
lsmask = combinestages(lsmask,"lsmask",Preferences,endsize,3);

% Calculate percentiles
p = [5,16,50,84,95];
gmst = cell2mat(cellfun(@(x) prctile(x, p), GMST, 'UniformOutput', false));

% Define eras
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
mesozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Selandian/Danian") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Induan"));
paleozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Tremadocian"));
paleonocopse = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));
allnocopse = find(GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));

% Calculate solar forcing
albedo_con = 0.29;
a = [0:ceil(GTS.Average(end))]';
b = linspace(0.29,0.26,numel(a));
albedo_lin = interp1(a,b,GTS.Average);
%albedo_lin = albedo_calc;
[CO2_corrected_con, dFtot_con, dFco2, dFsun_con, dFsa] = ...
    calcforcings(GTS.Average,PhanerozoicCO2,albedo_con);
[CO2_corrected_lin, dFtot_lin, ~, dFsun_lin, ~] = ...
    calcforcings(GTS.Average,PhanerozoicCO2,albedo_lin);
[CO2_corrected_all, dFtot_all, ~, dFsun_all, dFsa_all] = ...
    calcforcings(GTS.Average,PhanerozoicCO2,albedo_con,lsmask);

%% Plot Figure
figname = 'SupFig_Forcings.png';
fig = figure('Position',[-1290,41,665,700],'Color','w');
t = tiledlayout(3,5,'Padding','none','TileSpacing','compact');
cm = hex2rgb({'#004F60';'#0a9396';'#ffb703';'#ca6702';'#9b2226'},1);
% Panel 1: Forcings assuming constant solar
ax1 = nexttile([1,3]); hold on
plot(GTS.Average(allnocopse),median(dFco2(allnocopse,:),2),'Color',cm(2,:),'LineWidth',1.5);
plot(GTS.Average(allnocopse),median(dFsun_con(allnocopse,:),2),'Color',cm(3,:),'LineWidth',1.5);
plot(xlim,[0 0],'-','Color',cm(4,:),'LineWidth',1.5);
plot(GTS.Average(allnocopse),median(dFtot_con(allnocopse,:),2),'Color','k','LineWidth',2.5);
plot(xlim,[0 0],'k--')
ylim([-11,17])
geologictimescale(0,GTS.LowerBoundary(allnocopse(end)),'normal','reverse',gca,'standard','stages','off',10,2)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel('\DeltaF (Wm^{-2})','FontName','Arial','FontSize',13,'FontWeight','bold')
text(450,15.5,"A",'FontWeight','bold','FontName','Arial','FontSize',15,'Color','k')
text(450,13.5,"Constant α (0.29); no F\fontsize{9}geog",'FontWeight','bold','FontName','Arial','FontSize',13,'Color','k')
text(5,-3.5,"Total forcing",'FontWeight','bold','FontName','Arial','FontSize',11,'Color','k','HorizontalAlignment','right')
text(5,-5.5,"CO\fontsize{7}2 \fontsize{11}forcing",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(2,:),'HorizontalAlignment','right')
text(5,-7.5,"Solar forcing",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(3,:),'HorizontalAlignment','right')
text(5,-9.5,"Geographic forcing",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(4,:),'HorizontalAlignment','right')
% Panel 2: Constant crossplot
nexttile([1,2]); hold on, box on
X = prctile(dFtot_con,p,2);
Y = gmst;
xyerror(X(allnocopse,:),Y(allnocopse,:),p)
for ii = 1:size(X,1)
    if any(mesozoic == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#F8D99C',1),'MarkerSize',10);
    elseif any(cenozoic == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#66BBBD',1),'MarkerSize',10);
    elseif any(paleonocopse == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#DEA061',1),'MarkerSize',10);
    end
end
xlim([-13 18])
% ylim([6 45])
ylim([-10 45])
text(-12.25,42.5,"B",'FontWeight','bold','FontName','Arial','FontSize',15,'Color','k')
text(-12,40,sprintf("r = %.02f",corr(X(allnocopse,p==50),Y(allnocopse,p==50))),...
    'FontWeight','bold','FontName','Arial','FontSize',13,'Color','k')
text(8,14,"Cenozoic",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(2,:))
text(8,12.5,"Mesozoic",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(3,:))
text(8,11,"Paleozoic",'FontWeight','bold','FontName','Arial','FontSize',11,'Color',cm(4,:))
xlabel('\DeltaF (Wm^{-2})','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')

% Panel 3: Forcings assuming linear solar
nexttile([1,3]); hold on
plot(GTS.Average(allnocopse),median(dFco2(allnocopse,:),2),'Color',cm(2,:),'LineWidth',1.5);
plot(GTS.Average(allnocopse),median(dFsun_lin(allnocopse,:),2),'Color',cm(3,:),'LineWidth',1.5);
plot(GTS.Average(allnocopse),median(dFtot_lin(allnocopse,:),2),'Color','k','LineWidth',2.5);
plot(xlim,[0 0],'k--')
ylim([-11,17])
geologictimescale(0,GTS.LowerBoundary(allnocopse(end)),'normal','reverse',gca,'standard','stages','off',10,2)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel('\DeltaF (Wm^{-2})','FontName','Arial','FontSize',13,'FontWeight','bold')
text(450,15.5,"C",'FontWeight','bold','FontName','Arial','FontSize',15,'Color','k')
text(450,13.5,"Linear α (0.26 - 0.29); no F\fontsize{9}geog",'FontWeight','bold','FontName','Arial','FontSize',13,'Color','k')
% Panel 4: Linear crossplot
nexttile([1,2]); hold on, box on
X = prctile(dFtot_lin,p,2);
Y = gmst;
xyerror(X(allnocopse,:),Y(allnocopse,:),p)
for ii = 1:size(X,1)
    if any(mesozoic == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#F8D99C',1),'MarkerSize',10);
    elseif any(cenozoic == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#66BBBD',1),'MarkerSize',10);
    elseif any(paleonocopse == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#DEA061',1),'MarkerSize',10);
    end
end
xlim([-13 18])
%ylim([6 45])
ylim([-10 45])
text(-12.25,42.5,"D",'FontWeight','bold','FontName','Arial','FontSize',15,'Color','k')
text(-12,40,sprintf("r = %.02f",corr(X(allnocopse,p==50),Y(allnocopse,p==50))),...
    'FontWeight','bold','FontName','Arial','FontSize',13,'Color','k')
xlabel('\DeltaF (Wm^{-2})','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')

% Panel 5: Calcualted forcings
nexttile([1,3]); hold on
plot(GTS.Average(allnocopse),median(dFco2(allnocopse,:),2),'Color',cm(2,:),'LineWidth',1.5);
plot(GTS.Average(allnocopse),median(dFsun_all(allnocopse,:),2),'Color',cm(3,:),'LineWidth',1.5);
plot(GTS.Average(allnocopse),median(dFsa_all(allnocopse,:),2),'Color',cm(4,:),'LineWidth',1.5);
plot(GTS.Average(allnocopse),median(dFtot_all(allnocopse,:),2),'Color','k','LineWidth',2.5);
plot(xlim,[0 0],'k--')
ylim([-11,17])
geologictimescale(0,GTS.LowerBoundary(allnocopse(end)),'normal','reverse',gca,'standard','stages','off',10,2)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel('\DeltaF (Wm^{-2})','FontName','Arial','FontSize',13,'FontWeight','bold')
text(450,15.5,"E",'FontWeight','bold','FontName','Arial','FontSize',15,'Color','k')
text(450,13.5,"Constant α (0.29)",'FontWeight','bold','FontName','Arial','FontSize',13,'Color','k')
% Panel 6: Constant crossplot
nexttile([1,2]); hold on, box on
X = prctile(dFtot_all,p,2);
Y = gmst;
xyerror(X(allnocopse,:),Y(allnocopse,:),p)
for ii = 1:size(X,1)
    if any(mesozoic == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#F8D99C',1),'MarkerSize',10);
    elseif any(cenozoic == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#66BBBD',1),'MarkerSize',10);
    elseif any(paleonocopse == ii)
        plot(X(ii,p==50),Y(ii,p==50),'ko',...
            'MarkerFaceColor',hex2rgb('#DEA061',1),'MarkerSize',10);
    end
end
xlim([-13 18])
%ylim([6 45])
ylim([-10 45])
text(-12.25,42.5,"F",'FontWeight','bold','FontName','Arial','FontSize',15,'Color','k')
text(-12,40,sprintf("r = %.02f",corr(X(allnocopse,p==50),Y(allnocopse,p==50))),...
    'FontWeight','bold','FontName','Arial','FontSize',13,'Color','k')
xlabel('\DeltaF (Wm^{-2})','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')



if savefig
    export_fig(gcf,[figdir,'/Supplemental/',figname],'-p0.01','-m5')
end

%% Figure 2: Plot model albedo
figname = 'SupFig_Land.png';
fig = figure('Position',[-1500 330 730 316],'Color','w');
LS = lsstats(lsmask);
plot(GTS.Average,sum(LS.lssumweight,1)./sum(LS.lssumweight(:,1)),...
    '-','LineWidth',2,'Color',cm_col(3,:))
ylim([.4 1.25])
geologictimescale(0,GTS.LowerBoundary((end)),'normal','reverse',gca,'standard','stages','off',10,2)
ylabel('Continental area relative to today','FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')

if savefig
    export_fig(gcf,[figdir,'/Supplemental/',figname],'-p0.01','-m5')
end
