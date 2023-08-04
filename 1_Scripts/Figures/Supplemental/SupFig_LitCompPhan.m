%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Literature Comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
savefig = true;

% Select colormap
cm = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);

% PART 1: LOAD DATA
% Directory details
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
idx = contains(ItName,strcat("phCorr = ",string(pHCorr))) & ...
    contains(ItName,strcat("SnowballCorr = ",string(sbCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);

% Revise GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);

% SCOTESE
load("GMST_Scotese.mat")
bins = [GTS.UpperBoundary(1);GTS.LowerBoundary+.00000001];
binidx = discretize(GMST_Scotese.Age,bins);
GMSTscotese(:,1) = accumarray(binidx(~isnan(binidx)),GMST_Scotese.GMST(~isnan(binidx)),[],@mean);
GMSTscotese(GMSTscotese == 0) = NaN;

% LUNT
load("GMST_Lunt.mat")
x = 0:1:505;
y = interp1(GMST_Lunt.Time,GMST_Lunt.GMST,x);
binidx = discretize(x,bins)';
GMSTlunt(:,1) = accumarray(binidx(~isnan(binidx)),y(~isnan(binidx))',[],@mean);
GMSTlunt(GMSTlunt == 0) = NaN;

%% MAKE FIGURE
figname = 'SupFig_LitCompPhanerozoic.png';
fig = figure('Position',[-1248,491,479*2,420],'Color','w'); 
hold on, box on
GMSTmedian = cell2mat(cellfun(@(x) prctile(x,[5,16,50,84,95]), GMST, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[GMSTmedian(:,1);flipud(GMSTmedian(:,end))],'k','FaceColor',[.85 .85 .85],'EdgeColor','none');
fill([GTS.Average;flipud(GTS.Average)],[GMSTmedian(:,2);flipud(GMSTmedian(:,end-1))],'k','FaceColor',[.65 .65 .65],'EdgeColor','none');
p1 = plot(GTS.Average,GMSTmedian(:,3),'k-','LineWidth',2);
binidx = discretize(GMST_Scotese.Age,bins);
p2 = plot(GMST_Scotese.Age(~isnan(binidx)),GMST_Scotese.GMST(~isnan(binidx)),'-','LineWidth',1.5,'Color',cm(2,:));
p3 = plot(GTS.Average,GMSTscotese(:,1),'-','LineWidth',2,'Color',cm(1,:));
binidx = discretize(x,bins);
p4 = plot(x(~isnan(binidx)),y(~isnan(binidx)),'-','LineWidth',1.5,'Color',cm(end-1,:));
p5 = plot(GTS.Average,GMSTlunt(:,1),'-','LineWidth',2,'Color',cm(end,:));
pempty = plot(0,0,'-','Color','none');
geologictimescale(0,GTS.LowerBoundary(end),...
    'normal','reverse',gca,'standard','stages','off',9,2)
leg = legend([p1, pempty, p2, p3, p4, p5], {'PhanDA', '', ...
    'Scotese et al. (1 myr)','Scotese et al. (stage)',...
    'Valdes et al. (model timestep)','Valdes et al. (stage)'}, ...
    'NumColumns',3,'Location','north','FontName','Arial','FontSize',11);
set(gca,'FontName','Arial','FontSize',11)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')

% Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graveyard %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare: Scotese
Scotese = NaN(size(GTS,1),3);
for ii = 1:size(GMST,1)
    idx = find(GMST_Scotese.Age>GTS.UpperBoundary(ii) & ...
        GMST_Scotese.Age<=GTS.LowerBoundary(ii));
    if ~isempty(idx)
        Scotese(ii,:) = prctile(GMST_Scotese.GMST(idx),[5,50,95]);
    end
end
Scotese(1,:) = repmat(GMST_Scotese.GMST(1),1,3);
figurename = "LitComp_ScoteseLunt.png";
cm = hex2rgb({'#CA6702','#FFB703','#B4BE65','#E9D8A6','#9B2226','#0A9396'},1);
c=1; X = []; Y = [];
fig = figure('Position',[-1270,97,936,880]); 
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile
hold on, box on
for ii = 1:numel(Ref)
    if ~contains(Ref{ii}, "Anagnastou")
    for jj = 1:numel(Stage{ii})
        s = strsplit(Stage{ii}(jj),'/');
        idx = contains(GTS.Stage,s);
        ScotesePlot = mean(Scotese(idx,:),1);
        % Plot lit uncertainty
        plot(repmat(ScotesePlot(2),2,1),Temp{ii}{jj}(1,[1,3]),'-','color',[.5 .5 .5])
        % Plot PhanDA uncertainty
        plot(ScotesePlot([1,3]),repmat(Temp{ii}{jj}(1,2),2,1),'-','color',[.5 .5 .5])
        % Plot averages
        p{c} = plot(ScotesePlot(2),Temp{ii}{jj}(1,2),'ko','MarkerFaceColor',cm(ii,:));
        X(end+1,1) = ScotesePlot(2);
        Y(end+1,1) = Temp{ii}{jj}(1,2);
        if jj == 1
            c=c+1;
        end
        if contains(Ref{ii}, "Burls")
            % Plot lit uncertainty
            plot(repmat(ScotesePlot(2),2,1),Temp{ii}{jj}(2,[1,3]),'-','color',[.5 .5 .5])
            % Plot PhanDA uncertainty
            plot(ScotesePlot([1,3]),repmat(Temp{ii}{jj}(2,2),2,1),'-','color',[.5 .5 .5])
            % Plot averages
            p{c} = plot(ScotesePlot(2),Temp{ii}{jj}(2,2),'ks','MarkerFaceColor',cm(ii,:));
            X(end+1,1) = ScotesePlot(2);
            Y(end+1,1) = Temp{ii}{jj}(2,2);
            if jj == 1
                c=c+1;
            end
        end
    end   
    end
 end
mdl = fitlm(X,Y);
plot([10 40],[10,40],'-.','LineWidth',1.5,'Color',[.25,.25,.25])
leg = legend([p{5},p{1},p{2},p{3},p{4},p{6}], ...
    {Ref{5} + newline + "(95% CI)", ...
    Ref{1} + newline + "(90% CI)", ...
    Ref{2} + newline + "(Range)", ...
    Ref{3} + newline + "(all data; 2σ)", ...
    Ref{3} + newline + "(marine only; 2σ)", ...
    Ref{6}}, ...
    'Location','southeast','FontSize',12);
xlabel('Scotese et al. 2021 GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel('Time-slice GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
leg.FontSize = 11;

nexttile
hold on, box on
for ii = 1: size(d18bStage,1)-1
  % Plot lit uncertainty
    plot(repmat(Scotese(ii,2),2,1),d18bStage(ii,[1,3]),'-','color',[.5 .5 .5])
    % Plot PhanDA uncertainty
    plot(Scotese(ii,[1,3]),repmat(d18bStage(ii,2),2,1),'-','color',[.5 .5 .5])
end
plot(Scotese(1:ii,2),d18bStage(1:ii,2),'ko','MarkerFaceColor', ...
    hex2rgb('#005F73',1))
mdl = fitlm(Scotese(1:ii,2),d18bStage(1:ii,2));
x = [min(Scotese(1:ii,2)),max(Scotese(1:ii,2))]';
y = predict(mdl,x);
plot(x,y,'--','Color',hex2rgb('#005F73',1),'LineWidth',1.5)
set(gca,'Ydir','reverse')
xlabel('Scotese et al. 2021 GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel('\delta^{18}O_{benthic} (‰, VPDB)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylim([-2.8,5.7])
xlim([13.5 26])
set(fig,'Color','w')

%% Compare: Lunt
Lunt = NaN(size(GTS,1),3);
x = 0:1:65;
y = interp1(GMST_Lunt.Time,GMST_Lunt.GMST,x);
for ii = 1:size(GMST,1)
    idx = find(x>GTS.UpperBoundary(ii) & ...
        x<GTS.LowerBoundary(ii));
    if ~isempty(idx)
        Lunt(ii,:) = prctile(y(idx),[5,50,95]);
    end
end
Lunt(1,:) = repmat(y(1),1,3);
nexttile
hold on, box on
for ii = 1:numel(Ref)
    if ~contains(Ref{ii}, "Anagnastou")
    for jj = 1:numel(Stage{ii})
        s = strsplit(Stage{ii}(jj),'/');
        idx = contains(GTS.Stage,s);
        LuntPlot = mean(Lunt(idx,:),1);
        % Plot lit uncertainty
        plot(repmat(LuntPlot(2),2,1),Temp{ii}{jj}(1,[1,3]),'-','color',[.5 .5 .5])
        % Plot PhanDA uncertainty
        plot(LuntPlot([1,3]),repmat(Temp{ii}{jj}(1,2),2,1),'-','color',[.5 .5 .5])
        % Plot averages
        p{c} = plot(LuntPlot(2),Temp{ii}{jj}(1,2),'ko','MarkerFaceColor',cm(ii,:));
        X(end+1,1) = LuntPlot(2);
        Y(end+1,1) = Temp{ii}{jj}(1,2);
        if jj == 1
            c=c+1;
        end
        if contains(Ref{ii}, "Burls")
            % Plot lit uncertainty
            plot(repmat(LuntPlot(2),2,1),Temp{ii}{jj}(2,[1,3]),'-','color',[.5 .5 .5])
            % Plot PhanDA uncertainty
            plot(LuntPlot([1,3]),repmat(Temp{ii}{jj}(2,2),2,1),'-','color',[.5 .5 .5])
            % Plot averages
            p{c} = plot(LuntPlot(2),Temp{ii}{jj}(2,2),'ks','MarkerFaceColor',cm(ii,:));
            X(end+1,1) = LuntPlot(2);
            Y(end+1,1) = Temp{ii}{jj}(2,2);
            if jj == 1
                c=c+1;
            end
        end
    end   
    end
 end
mdl = fitlm(X,Y);
% Y = predict(mdl,X);
% plot(X,Y,'k--','LineWidth',1.5)
plot([12 32],[12 32],'-.','LineWidth',1.5,'Color',[.25,.25,.25])
xlabel('Model-only GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel('Time-slice GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
leg.FontSize = 11;


nexttile
hold on, box on
for ii = 1: size(d18bStage,1)
  % Plot lit uncertainty
    plot(repmat(Lunt(ii,2),2,1),d18bStage(ii,[1,3]),'-','color',[.5 .5 .5])
    % Plot PhanDA uncertainty
    plot(Lunt(ii,[1,3]),repmat(d18bStage(ii,2),2,1),'-','color',[.5 .5 .5])
end
plot(Lunt(1:ii,2),d18bStage(1:ii,2),'ko','MarkerFaceColor', ...
    hex2rgb('#005F73',1))
mdl = fitlm(Lunt(1:ii,2),d18bStage(1:ii,2));
x = [min(Lunt(1:ii,2)),max(Lunt(1:ii,2))]';
y = predict(mdl,x);
plot(x,y,'--','Color',hex2rgb('#005F73',1),'LineWidth',1.5)
set(gca,'Ydir','reverse')
xlabel('Model-only GMST (^oC)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylabel('\delta^{18}O_{benthic} (‰, VPDB)','FontName','Arial','FontSize',15,'FontWeight','bold')
ylim([-2.8,5.7])
xlim([13.5 26])
set(fig,'Color','w')
export_fig(fig,figurename,'-p0.01','-m5')

