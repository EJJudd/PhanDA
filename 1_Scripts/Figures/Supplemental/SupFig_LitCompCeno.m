%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Cenozoic Literature Comparison %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
savefig = true;

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")

% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
swCorr = ["snowball", "off"];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",swCorr)) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);
% Revise GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
P = [5 16 50 84 95];
pgmst = cell2mat(cellfun(@(x) prctile(x,P), GMST, 'UniformOutput', false));
[Ref, Stage, Age, Temp, PhanDA, Lit, Ridx] = loadcenolitdata(GMST,GTS);

% Benthic stack
%Load data
d18b = readtable('/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Data/SupplementalData/d18Obenthic_Westerhold.csv');
d18b(d18b.Time>66,:) = [];
load("GMST_Hansen.mat");
load("PhanerozoicpHv6.mat");
% Assign age bins
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
bins = [GTS.UpperBoundary(1);GTS.LowerBoundary(cenozoic)+.00000001];
P = [5 16 50 84 95];
d18bStage = NaN(22,numel(P));
GMSThansen = NaN(22,numel(P));
GMSThansenpH = NaN(22,numel(P));

% Calculate the pH corrected value
pHbin = discretize(d18b.Time,bins);
d18b.pH = median(PhanerozoicpH(pHbin,:),2);
d18b.d18Oph = d18b.ISOBENd18oLOESSsmooth - (1.42* (d18b.pH(1) - d18b.pH));

% define the holocene, lgm & oi-1 d18Ob
hol = 3.88;
lgm = 5.3;
oi1 = 1.5;
% % Calculate bottom water temp
d18b.Tb = zeros(size(d18b.Time));
% if d18Ob <= 1.5
d18b.Tb(d18b.ISOBENd18oLOESSsmooth<=1.5) = ...
    -4*d18b.d18Oph(d18b.ISOBENd18oLOESSsmooth<=1.5) + 12;
% if d18Ob > 3.88
d18b.Tb(d18b.ISOBENd18oLOESSsmooth>3.88) = ...
    1-1.41*(d18b.d18Oph(d18b.ISOBENd18oLOESSsmooth>3.88) - 3.88);
% if d18Ob >1.5 * <= 3.88
d18b.Tb(d18b.Tb == 0) = ...
    6-2.1*(d18b.d18Oph(d18b.Tb == 0) - 1.5);
% Calcate GMST
d18b.GMST = d18b.Tb - 0.35 * (exp(0.8*(d18b.d18Oph - oi1)) - 1) + 15;


% by Stage
for ii = 1:numel(cenozoic)
    idx = find(d18b.Time>=GTS.UpperBoundary(ii) & ...
        d18b.Time<GTS.LowerBoundary(ii));
    if ~isempty(idx)
        d18bStage(ii,1:5) = prctile(d18b.ISOBENd18oLOESSsmooth(idx),P);
        GMSThansen(ii,1:5) = prctile(GMST_Hansen.GMST(idx),P);
        GMSThansenpH(ii,1:5) = prctile(d18b.GMST(idx),P);
    end   
end




%% MAKE FIGURES

%% FIGURE 1: LIT COMP
figname = 'SupFig_LitCompCenozoic.png';
figure('Position',[-1636 550 946 332],'Color','w')
tiledlayout(1,9,'Padding','none','TileSpacing','compact');
% fig = figure('Position',[-1248,491,650,420],'Color','w'); hold on, box on
%cm = hex2rgb({'#CA6702','#FFB703','#0A9396','#E9D8A6','#9B2226','#005f73','#B4BE65'},1);
cm = hex2rgb({'#005f73','#0a9396','#94d2bd','#b4be65','#ffb703','#ca6702','#9b2226'},1);
[reflist,nr] = unique(Ridx,'stable');
reflist = reflist([6,7,3,4,2,8,1,5]);
reforder = string(Ref([6,7,3,2,8,1,5]));
ax1 = nexttile([1,5]); hold on
GMSTmedian = cell2mat(cellfun(@(x) prctile(x,[5,16,50,84,95]), GMST, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[GMSTmedian(:,1);...
    flipud(GMSTmedian(:,end))],'k','FaceColor',[.85 .85 .85],'EdgeColor','none');
fill([GTS.Average;flipud(GTS.Average)],...
    [GMSTmedian(:,2);flipud(GMSTmedian(:,end-1))],'k','FaceColor',...
    [.65 .65 .65],'EdgeColor','none');
plot(GTS.Average,GMSTmedian(:,3),'k-','LineWidth',2);
ylim([10 43])
geologictimescale(0,GTS.LowerBoundary(22),'normal',...
    'reverse',ax1,'standard','stages','off',8,1)
xlabel('Age (Ma)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel(['Tropical SST (',char(176),'C)'],'FontName','Arial','FontSize',15,'FontWeight','bold')


ax2 = nexttile([1,4]); hold on, box on
for ii = 1:numel(reflist)
    % Find data
    idx = find(Ridx == reflist(ii));
    % Assign plot color
    pc = find(reforder == strrep(strrep(reflist(ii),' all',''),' sst',''));
    sa = find(string(Ref) == strrep(strrep(reflist(ii),' all',''),' sst',''));
    % Plot uncertainty
    errorbar(ax2,PhanDA(idx,2),Lit(idx,2),Lit(idx,2)-Lit(idx,1),...
        Lit(idx,3)-Lit(idx,2),PhanDA(idx,2)-PhanDA(idx,1),...
        PhanDA(idx,3)-PhanDA(idx,2),'.','Color',[.5,.5,.5],...
        'LineWidth',1,'CapSize',0)
    % Plot age uncertainty
    x = NaN(size(Age{sa})); y = NaN(size(Age{sa}));
    % Plot averages
    if ~contains(reflist(ii),'sst')
        for jj = 1:size(Age{sa},1)
            x(jj) = mean(Age{sa}{jj});
            y(jj) = Temp{sa}{jj}(1,2);
            plot(ax1,[min(Age{sa}{jj}),max(Age{sa}{jj})],repmat(y(jj),1,2),'k-')
            plot(ax1,repmat(x(jj),1,2),[Temp{sa}{jj}(1,1),Temp{sa}{jj}(1,3)],'k-')
        end
            p{ii} = plot(ax2,PhanDA(idx,2),Lit(idx,2),'ko','MarkerFaceColor',cm(pc,:),...
            'MarkerSize',11);
        plot(ax1,x,y,'ko','MarkerFaceColor',cm(pc,:),'MarkerSize',12)
    else
        for jj = 1:size(Age{sa},1)
            x(jj) = mean(Age{sa}{jj});
            y(jj) = Temp{sa}{jj}(2,2);
            plot(ax1,[min(Age{sa}{jj}),max(Age{sa}{jj})],repmat(y(jj),1,2),'k-')
            plot(ax1,repmat(x(jj),1,2),[Temp{sa}{jj}(2,1),Temp{sa}{jj}(2,3)],'k-')
        end
        p{ii} = plot(ax2,PhanDA(idx,2),Lit(idx,2),'ks','MarkerFaceColor',...
            cm(pc,:),'MarkerSize',11);
        plot(ax1,x,y,'ks','MarkerFaceColor',cm(pc,:),'MarkerSize',13)
    end
end

% Record stages
Stages = Stage([1:3,5:end]);
stages = [];
for ii = 1:numel(Stages)
    stages = [stages;Stages{ii}];
end
stages = unique(stages);

ylim([12 43]), xlim(ylim)
plot(ylim,ylim,'--','LineWidth',1.5,'Color','k')
leg = legend([p{end},p{end-1},p{end-2},p{end-3},p{end-5},p{end-4},p{end-6},p{end-7}], ...
    {reforder(end) + newline + "(95% CI)", ...
    reforder(end-1) + newline + "(90% CI)", ...
    reforder(end-2) + newline + "(95% CI)", ...
    reforder(end-3) + newline + "(Range)", ...
    reforder(end-4) + newline + "(all data; 2σ)", ...
    reforder(end-4) + newline + "(marine only; 2σ)", ...
    reforder(end-5) + newline + "(Range)", ...
    reforder(end-6) + newline + "(95% CI)"}, ...
    'Location','eastoutside','FontSize',12,'NumColumns',1);
ax = gca;
ax.FontSize = 11;
xlabel(['PhanDA GMST (',char(176),'C)'],'FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel(['Literature GMST (',char(176),'C)'],'FontName','Arial','FontSize',15,'FontWeight','bold')
leg.FontSize = 11;
text(13,40.5,sprintf('r = %.2f',corr(PhanDA(:,2),Lit(:,2))),...
    'FontName','Arial','FontSize',13,'FontWeight','bold')

text(ax1, 65,41.5,'A','FontName','Arial','FontSize',15,'FontWeight','bold')
text(ax2, 13,41.5,'B','FontName','Arial','FontSize',15,'FontWeight','bold')

% Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end

%% FIGURE 2: BENTHIC COMP
figname = 'SupFig_BenthicStack.png';
cm = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
gmst = cell2mat(cellfun(@(x) prctile(x,P), GMST, 'UniformOutput', false));

fig = figure('Units','inches','Position',[-22,6,8,7.25],'Color','w');
tiledlayout(2,2,'Padding','none','TileSpacing','compact');

ax1 = nexttile(1); hold on,
fill([GTS.Average(cenozoic);flipud(GTS.Average(cenozoic))],[gmst(cenozoic,1);flipud(gmst(cenozoic,end))],'k','FaceColor',[.85 .85 .85],'EdgeColor','none');
fill([GTS.Average(cenozoic);flipud(GTS.Average(cenozoic))],[gmst(cenozoic,2);flipud(gmst(cenozoic,end-1))],'k','FaceColor',[.65 .65 .65],'EdgeColor','none');
p1 = plot(GTS.Average(cenozoic),gmst(cenozoic,3),'k-','LineWidth',2);
ax1.TickLength = [.0175 .025];
ylim([0 40]),set(ax1,'YTick',[10:10:40])
geologictimescale(0,GTS.LowerBoundary(cenozoic(end)),...
    'normal','reverse',gca,'standard','stages','off',7,2)
yyaxis right
p2 = plot(d18b.Time,d18b.ISOBENd18oLOESSsmooth,'-','Color',cm(2,:));
p3 = plot(GTS.Average(cenozoic),d18bStage(cenozoic,3),'-','LineWidth',2,'color',cm(1,:));
ylim(ax1,[-6 8])
set(ax1,'YDir','reverse','YTick',[-2:2:6])
ax1.YColor = cm(2,:);
ylabel(ax1,['\fontsize{13}δ^{18}O','\fontsize{7}benthic', '\fontsize{13}  (‰)           '],'FontName','Arial','FontSize',13,'FontWeight','bold')
yyaxis left
ylabel(ax1,['                 GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
text(ax1,64.5,36,'A','FontName','Arial','FontSize',15,'FontWeight','bold')

ax2 = nexttile;hold on, box on
for ii = 1:cenozoic(end)
  % Plot lit uncertainty
    plot(repmat(gmst(ii,3),2,1),d18bStage(ii,[1,end]),'-','color',[.5 .5 .5])
    % Plot PhanDA uncertainty
    plot([gmst(ii,1),gmst(ii,end)],repmat(d18bStage(ii,3),2,1),'-','color',[.5 .5 .5])
end
scatter(gmst(cenozoic,3),d18bStage(:,3),75,GTS.Average(cenozoic),'filled', ...
    'MarkerEdgeColor','k')
c = colorbar('south');
caxis([0, GTS.Average(cenozoic(end))])
colormap(viridis)
% York
sigma_X = (gmst(cenozoic,4) - gmst(cenozoic,2))./2;
sigma_Y = (d18bStage(:,4) - d18bStage(:,2))./2;
r = 0;
[b, m, sb, sm] = york_fit(gmst(cenozoic,3)',d18bStage(:,3)',sigma_X',sigma_Y', r);
x = [min(gmst(:,3)),max(gmst(:,3))];
y = m.*x + b;
plot(x,y,'k-','LineWidth',1.5)
set(gca,'Ydir','reverse')
xlabel(c,'Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\fontsize{13}δ^{18}O','\fontsize{7}benthic', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylim([-2.5,5.9])
xlim([9,40])
text(ax2,9.75,-2.1,'B','FontName','Arial','FontSize',15,'FontWeight','bold')
text(ax2,9.75,-1.7,sprintf('r = %.2f',corr(gmst(cenozoic,3),d18bStage(:,3))),...
    'FontName','Arial','FontSize',13,'FontWeight','bold')

ax3 = nexttile;hold on, box on
n = size(GMSThansen,1);
errorbar(pgmst(1:n,3),GMSThansen(:,3),...
    GMSThansen(:,3)-GMSThansen(:,1),...
    GMSThansen(:,5)-GMSThansen(:,3),...
    pgmst(1:n,3)-pgmst(1:n,1),...
    pgmst(1:n,5)-pgmst(1:n,3),'.','Color',[.5,.5,.5],...
        'LineWidth',1,'CapSize',0)
scatter(pgmst(1:n,3),GMSThansen(:,3),75,GTS.Average(cenozoic),'filled', ...
    'MarkerEdgeColor','k')
xlim(ylim)
ax3.YTick = [10:5:40];
plot(xlim,xlim,'k--','LineWidth',1)
sigma_X = (gmst(cenozoic,4) - gmst(cenozoic,2))./2;
sigma_Y = (GMSThansen(:,4) - GMSThansen(:,2))./2;
r = 0;
[b, m, sb, sm] = york_fit(gmst(cenozoic,3)',GMSThansen(:,3)',sigma_X',sigma_Y', r);
x = [min(gmst(:,3)),max(gmst(:,3))];
y = m.*x + b;
plot(x,y,'k','LineWidth',1.5)
xlabel(['PhanDA GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['Scaled GMST',newline,'(',char(176),'C; Hansen et al. 2023)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
text(ax3,5.5,38.5,'C','FontName','Arial','FontSize',15,'FontWeight','bold')
text(ax3,27.5,20,sprintf('slope = %.1f',m),...
    'FontName','Arial','FontSize',13,'FontWeight','bold')


ax4 = nexttile;hold on, box on
n = size(GMSThansenpH,1);
errorbar(pgmst(1:n,3),GMSThansenpH(:,3),...
    GMSThansenpH(:,3)-GMSThansenpH(:,1),...
    GMSThansenpH(:,5)-GMSThansenpH(:,3),...
    pgmst(1:n,3)-pgmst(1:n,1),...
    pgmst(1:n,5)-pgmst(1:n,3),'.','Color',[.5,.5,.5],...
        'LineWidth',1,'CapSize',0)
scatter(pgmst(1:n,3),GMSThansenpH(:,3),75,GTS.Average(cenozoic),'filled', ...
    'MarkerEdgeColor','k')
xlim([5,40])
ylim(xlim)
plot(xlim,xlim,'k--','LineWidth',1)
ax3.YTick = [10:5:40];
sigma_X = (gmst(cenozoic,4) - gmst(cenozoic,2))./2;
sigma_Y = (GMSThansenpH(:,4) - GMSThansenpH(:,2))./2;
r = 0;
[b, m, sb, sm] = york_fit(gmst(cenozoic,3)',GMSThansenpH(:,3)',sigma_X',sigma_Y', r);
x = [min(gmst(:,3)),max(gmst(:,3))];
y = m.*x + b;
plot(x,y,'k','LineWidth',1.5)
xlabel(['PhanDA GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['pH corrected scaled GMST',newline,'(',char(176),'C; Hansen et al. 2023)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
text(ax4,5.5,38.5,'D','FontName','Arial','FontSize',15,'FontWeight','bold')
text(ax4,27.5,20,sprintf('slope = %.1f',m),...
    'FontName','Arial','FontSize',13,'FontWeight','bold')

fprintf('r = %.2f\n',corr(gmst(cenozoic,3),d18bStage(:,3)))
fprintf('r = %.2f\n',corr(gmst(cenozoic,3),GMSThansenpH(:,3)))
fprintf('r = %.2f\n',corr(gmst(cenozoic,3),GMSThansen(:,3)))

% Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end
