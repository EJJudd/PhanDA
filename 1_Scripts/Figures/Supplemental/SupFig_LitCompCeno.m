%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Cenozoic Literature Comparison %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
savefig = true;

% PART 1: LOAD DATA
% Directory details
% assdate = '03Jul2023';
assdate = '27Jul2023';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","ItName")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")

% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
sbCorr = [true, false];
rMeth = ["low","medium","high"];
%rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",string(pHCorr))) & ...
    contains(ItName,strcat("SnowballCorr = ",string(sbCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);
% Revise GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);

[Ref, Stage, Age, Temp, PhanDA, Lit, Ridx] = loadcenolitdata(GMST,GTS);


% Benthic stack
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
bins = [GTS.UpperBoundary(1);GTS.LowerBoundary(cenozoic)+.00000001];
d18b = readtable('/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Data/SupplementalData/d18Obenthic_Westerhold.csv');
P = [5 16 50 84 95];
d18bStage = NaN(22,numel(P));
% by Stage
for ii = 1:numel(cenozoic)
    idx = find(d18b.Time>=GTS.UpperBoundary(ii) & ...
        d18b.Time<GTS.LowerBoundary(ii));
    if ~isempty(idx)
        d18bStage(ii,1:5) = prctile(d18b.ISOBENd18oLOESSsmooth(idx),P);
    end   
end


%% MAKE FIGURES

%% FIGURE 1: LIT COMP
figname = 'SupFig_LitCompCenozoic.png';
fig = figure('Position',[-1248,491,650,420],'Color','w'); hold on, box on
%cm = hex2rgb({'#CA6702','#FFB703','#0A9396','#E9D8A6','#9B2226','#005f73','#B4BE65'},1);
cm = hex2rgb({'#005f73','#0a9396','#94d2bd','#b4be65','#ffb703','#ca6702','#9b2226'},1);
[reflist,nr] = unique(Ridx,'stable');
reflist = reflist([6,7,3,4,2,8,1,5]);
reforder = string(Ref([6,7,3,2,8,1,5]));
c = 0;
for ii = 1:numel(reflist)
    % Find data
    idx = find(Ridx == reflist(ii));
    % Assign plot color
    pc = find(reforder == strrep(strrep(reflist(ii),' all',''),' sst',''));
    % Plot uncertainty
    errorbar(PhanDA(idx,2),Lit(idx,2),Lit(idx,2)-Lit(idx,1),...
        Lit(idx,3)-Lit(idx,2),PhanDA(idx,2)-PhanDA(idx,1),...
        PhanDA(idx,3)-PhanDA(idx,2),'.','Color',[.5,.5,.5],...
        'LineWidth',1,'CapSize',0)
    % Plot averages
    if ~contains(reflist(ii),'sst')
        p{ii} = plot(PhanDA(idx,2),Lit(idx,2),'ko','MarkerFaceColor',cm(pc,:),...
            'MarkerSize',11);
    else
        p{ii} = plot(PhanDA(idx,2),Lit(idx,2),'ks','MarkerFaceColor',...
            cm(pc,:),'MarkerSize',11);
    end
end

ylim([12 43]), xlim(ylim)
plot(ylim,ylim,'--','LineWidth',1.5,'Color','k')
leg = legend([p{end},p{end-1},p{end-2},p{end-3},p{end-4},p{end-5},p{end-6},p{end-7}], ...
    {reforder(end) + newline + "(95% CI)", ...
    reforder(end-1) + newline + "(90% CI)", ...
    reforder(end-2) + newline + "(95% CI)", ...
    reforder(end-3) + newline + "(Range)", ...
    reforder(end-4) + newline + "(Range)", ...
    reforder(end-4) + newline + "(all data; 2σ)", ...
    reforder(end-5) + newline + "(marine only; 2σ)", ...
    reforder(end-6)}, ...
    'Location','eastoutside','FontSize',12,'NumColumns',1);
ax = gca;
ax.FontSize = 11;
xlabel(['PhanDA GMST (',char(176),'C)'],'FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel(['Literature GMST (',char(176),'C)'],'FontName','Arial','FontSize',15,'FontWeight','bold')
leg.FontSize = 11;
text(13,41,sprintf('r = %.2f',corr(PhanDA(:,2),Lit(:,2))),...
    'FontName','Arial','FontSize',13,'FontWeight','bold')

% Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end

%% FIGURE 2: BENTHIC COMP
figname = 'SupFig_BenthicStack.png';
cm = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
gmst = cell2mat(cellfun(@(x) prctile(x,P), GMST, 'UniformOutput', false));

fig = figure('Units','inches','Position',[-22,6,10,5],'Color','w');
tiledlayout(2,2,'Padding','none','TileSpacing','compact');

ax1 = nexttile(1); hold on, box on
p2 = plot(d18b.Time,d18b.ISOBENd18oLOESSsmooth,'-','Color',cm(2,:));
p3 = plot(GTS.Average(cenozoic),d18bStage(cenozoic,3),'-','LineWidth',2,'color',cm(1,:));
ylim(ax1,[-3 5.5])
set(ax1,'YDir','reverse','XDir','reverse','xlim',[0,GTS.LowerBoundary(cenozoic(end))],'XMinorTick','on')
ax1.TickLength = [.0175 .025];
ylabel(['\fontsize{13}δ^{18}O','\fontsize{7}benthic', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
text(ax1,64.5,-2.25,'A','FontName','Arial','FontSize',15,'FontWeight','bold')

ax2 = nexttile(3); hold on,
fill([GTS.Average(cenozoic);flipud(GTS.Average(cenozoic))],[gmst(cenozoic,1);flipud(gmst(cenozoic,end))],'k','FaceColor',[.85 .85 .85],'EdgeColor','none');
fill([GTS.Average(cenozoic);flipud(GTS.Average(cenozoic))],[gmst(cenozoic,2);flipud(gmst(cenozoic,end-1))],'k','FaceColor',[.65 .65 .65],'EdgeColor','none');
p1 = plot(GTS.Average(cenozoic),gmst(cenozoic,3),'k-','LineWidth',2);
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ax2.TickLength = [.0175 .025];
ylim(ax2,[10 41])
geologictimescale(0,GTS.LowerBoundary(cenozoic(end)),...
    'normal','reverse',gca,'standard','stages','off',7,2)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold'), box on
text(ax2,64.5,38,'B','FontName','Arial','FontSize',15,'FontWeight','bold')
ax3 = nexttile([2,1]);hold on, box on
for ii = 1:cenozoic(end)
  % Plot lit uncertainty
    plot(repmat(gmst(ii,3),2,1),d18bStage(ii,[2,end-1]),'-','color',[.5 .5 .5])
    % Plot PhanDA uncertainty
    plot([gmst(ii,2),gmst(ii,end-1)],repmat(d18bStage(ii,3),2,1),'-','color',[.5 .5 .5])
end
scatter(gmst(cenozoic,3),d18bStage(:,3),75,GTS.Average(cenozoic),'filled', ...
    'MarkerEdgeColor','k')
c = colorbar('south');
caxis([0, GTS.Average(cenozoic(end))])
% York
sigma_X = (gmst(cenozoic,4) - gmst(cenozoic,2))./2;
sigma_Y = (d18bStage(:,4) - d18bStage(:,2))./2;
r = 0;
[b, m, sb, sm] = york_fit(gmst(cenozoic,3)',d18bStage(:,3)',sigma_X',sigma_Y', r);
x = [min(gmst(:,3)),max(gmst(:,3))];
y = m.*x + b;
plot(x,y,'k--','LineWidth',1.5)
set(gca,'Ydir','reverse')
xlabel(c,'Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\fontsize{13}δ^{18}O','\fontsize{7}benthic', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylim([-2,5.9])
xlim([9,37])
text(ax3,9.5,-1.6,'C','FontName','Arial','FontSize',15,'FontWeight','bold')
text(ax3,9.5,-1.15,sprintf('r = %.2f',corr(gmst(cenozoic,3),d18bStage(:,3))),...
    'FontName','Arial','FontSize',13,'FontWeight','bold')

fprintf('r = %.2f\n',corr(gmst(cenozoic,3),d18bStage(:,3)))

% Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graveyard %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BENTHIC STACK
d18b = readtable('/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Data/SupplementalData/d18Obenthic_Westerhold.csv');
% by Stage
for ii = 1:22
    idx = find(d18b.Time>=GTS.UpperBoundary(ii) & ...
        d18b.Time<GTS.LowerBoundary(ii));
    if ~isempty(idx)
        d18bStage(ii,1:3) = prctile(d18b.ISOBENd18oLOESSsmooth(idx),[5 50 95]);
        d18bStage(ii,4) = std(d18b.ISOBENd18oLOESSsmooth(idx));
    end   
end

nexttile
hold on, box on
for ii = 1: size(d18bStage,1)
  % Plot lit uncertainty
    plot(repmat(median(GMST(ii,:)),2,1),d18bStage(ii,[1,3]),'-','color',[.5 .5 .5])
    % Plot PhanDA uncertainty
    plot(prctile(GMST(ii,:),[5,95]),repmat(d18bStage(ii,2),2,1),'-','color',[.5 .5 .5])
end
plot(median(GMST(1:ii,:),2),d18bStage(:,2),'ko','MarkerFaceColor', ...
    hex2rgb('#005F73',1))
mdl = fitlm(median(GMST(1:ii,:),2),d18bStage(:,2));
x = [min(median(GMST(1:ii,:),2)),max(median(GMST(1:ii,:),2))]';
y = predict(mdl,x);
plot(x,y,'--','Color',hex2rgb('#005F73',1),'LineWidth',1.5)
set(gca,'Ydir','reverse')
xlabel('PhanDA GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel('\delta^{18}O_{benthic} (‰, VPDB)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylim([-2.8,5.7])
xlim([9,39])
set(fig,'Color','w')
export_fig(fig,figname,'-p0.01','-m5')


%% Compare - PhanDA on Y axis
cm = hex2rgb({'#CA6702','#FFB703','#B4BE65','#E9D8A6','#9B2226','#0A9396'},1);
c=1;
fig = figure('Position',[-1248,491,936,420]); 
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile
hold on, box on
for ii = 1:numel(Ref)
    if ~contains(Ref{ii}, "Anagnastou")
    for jj = 1:numel(Stage{ii})
        s = strsplit(Stage{ii}(jj),'/');
        idx = contains(GTS.Stage,s);
        PhanDA = prctile(GMST(idx,:),[5,50,95],'all');
        % Plot lit uncertainty
        plot(Temp{ii}{jj}(1,[1,3]),repmat(PhanDA(2),2,1),'-','color',[.5 .5 .5])
        % Plot PhanDA uncertainty
        plot(repmat(Temp{ii}{jj}(1,2),2,1),PhanDA([1,3]),'-','color',[.5 .5 .5])
        % Plot averages
        p{c} = plot(Temp{ii}{jj}(1,2),PhanDA(2),'ko','MarkerFaceColor',cm(ii,:));
        if jj == 1
            c=c+1;
        end
        if contains(Ref{ii}, "Burls")
            % Plot lit uncertainty
            plot(Temp{ii}{jj}(2,[1,3]),repmat(PhanDA(2),2,1),'-','color',[.5 .5 .5])
            % Plot PhanDA uncertainty
            plot(repmat(Temp{ii}{jj}(2,2),2,1),PhanDA([1,3]),'-','color',[.5 .5 .5])
            % Plot averages
            p{c} = plot(Temp{ii}{jj}(2,2),PhanDA(2),'ks','MarkerFaceColor',cm(ii,:));
            if jj == 1
                c=c+1;
            end
        end
    end   
    end
end
plot([10 40],[10,40],'k-.','LineWidth',1.5)
leg = legend([p{5},p{1},p{2},p{3},p{4},p{6}], ...
    {Ref{5} + newline + "(95% CI)", ...
    Ref{1} + newline + "(90% CI)", ...
    Ref{2} + newline + "(Range)", ...
    Ref{3} + newline + "(all data; 2σ)", ...
    Ref{3} + newline + "(marine only; 2σ)", ...
    Ref{6}}, ...
    'Location','northwest','FontSize',12);
ylabel('PhanDA GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
xlabel('Literature GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
leg.FontSize = 11;

nexttile
hold on, box on
for ii = 1: size(d18bStage,1)
  % Plot lit uncertainty
    plot(d18bStage(ii,[1,3]),repmat(median(GMST(ii,:)),2,1),'-','color',[.5 .5 .5])
    % Plot PhanDA uncertainty
    plot(repmat(d18bStage(ii,2),2,1),prctile(GMST(ii,:),[5,95]),'-','color',[.5 .5 .5])
end
plot(d18bStage(:,2),median(GMST(1:ii,:),2),'ko','MarkerFaceColor', ...
    hex2rgb('#005F73',1))
% mdl = fitlm(d18bStage(:,2),median(GMST(1:ii,:),2));
% x = [min(d18bStage(:,2))-.25,max(d18bStage(:,2))+.25]';
% y = predict(mdl,x);
% plot(x,y,'--','Color',hex2rgb('#005F73',1),'LineWidth',1.5)
set(gca,'Xdir','reverse')
ylabel('PhanDA GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
xlabel('\delta^{18}O_{benthic} (‰, VPDB)','FontName','Arial','FontSize',15,'FontWeight','bold')
xlim([-2.8,5.7])
ylim([9,39])
set(fig,'Color','w')

k = normrnd(repmat(d18bStage(:,2),1,1000),repmat(d18bStage(:,4),1,1000));
y = itfit(k,GMST(1:22,:),1000,[min(d18bStage(:,1))+1;d18bStage(:,2);max(d18bStage(:,3))-.25]);
[xplot,idx] = sort([min(d18bStage(:,1))+1;d18bStage(:,2);max(d18bStage(:,3))-.25]);
yplot = y(idx,:);
f = fill([xplot;flipud(xplot)],[prctile(yplot,5,2);flipud(prctile(yplot,95,2))],'k');
f.FaceAlpha = .3; f.EdgeAlpha = 0; f.FaceColor = hex2rgb('#005F73',1);
plot(xplot,median(yplot,2),'--','LineWidth',2,'Color',hex2rgb('#005F73',1))

