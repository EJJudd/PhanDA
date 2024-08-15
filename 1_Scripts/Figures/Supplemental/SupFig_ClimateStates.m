%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Climate states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
savefig = true;

% Select colormap
cm = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","LTG","ItName","TASprior")
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
LTG = combineruns(LTG,idx,2);

% Revise LTG, GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
LTG = combinestages(LTG,"LTG",Preferences,endsize,2);

% Calculate prior GMST and LTG AND remove duplicates
GMSTprior = cell2mat(latweightgmst(TASprior));
LTGprior = cellfun(@(x) squeeze(mean(x, 2)), TASprior, 'UniformOutput', false);
LTGprior = horzcat(LTGprior{:});
u = unique([GMSTprior,LTGprior'],'rows','stable');
GMSTprior = u(:,1);
LTGprior = u(:,2:end)';

% Define the regions
poles = abs(Lat)>=66.5;
tropics = abs(Lat)<=23.5;

% subdivide GMST of **posterior** percentile
prctiles = [20:20:100];
mgmst = cell2mat(cellfun(@(x) nanmedian(x), GMST, 'UniformOutput', false));
Ppost = prctile(mgmst, prctiles);
Pprior = prctile(GMSTprior, prctiles);
%Pprior = prctile(GMSTprior, prctiles);
CSpoidx.ih = find(mgmst<=Ppost(1));
CSpoidx.ch = find(mgmst>Ppost(1) & mgmst<=Ppost(2));
CSpoidx.tr = find(mgmst>Ppost(2) & mgmst<=Ppost(3));
CSpoidx.gh = find(mgmst>Ppost(3) & mgmst<=Ppost(4));
CSpoidx.hh = find(mgmst>Ppost(4));


%% PART 3: MAKE FIGURES
% FIGURE 1 - summary prior LTG
fig = figure('Units','inches','Position',[-22,6,12.875,4.2],'Color','w');
tiledlayout(1,3,'Padding','none','TileSpacing','compact');
% Panel A - Posterior GMST vs tropical-polar range
ax1 = nexttile;
hold on, box on
fn = fieldnames(CSpoidx);
X = NaN(numel(GMST),3);
Y = NaN(numel(GMST),3);
Yt = NaN(numel(GMST),3);
Yp = NaN(numel(GMST),3);
meandim = 1;
for ii = 1:numel(GMST)
    yt = LTG{ii}; yt(~tropics,:) = NaN;
    yp = LTG{ii}; yp(~poles,:) = NaN;
    X(ii,:) = prctile(GMST{ii},[16,50,84]);
    Y(ii,:) = prctile(latweightgmst(yt,meandim)-latweightgmst(yp,meandim),...
        [16,50,84]);
    Yt(ii,:) = prctile(latweightgmst(yt,meandim),[16,50,84]);
    Yp(ii,:) = prctile(latweightgmst(yp,meandim),[16,50,84]);
end
%Y = GMSTgrad(:,2:4);
errorbar(X(:,2),Y(:,2),Y(:,2)-Y(:,1),Y(:,3)-Y(:,2),X(:,2)-X(:,1),X(:,3)-X(:,2),...
    '.','Color',[.75 .75 .75],'CapSize',0)
for ii = 1:numel(fn)
    plot(X(CSpoidx.(fn{ii}),2),Y(CSpoidx.(fn{ii}),2),'ko','MarkerFaceColor',cm(ii,:),'MarkerSize',10)
end
ax1.FontName = 'Arial';ax1.FontSize = 11;
xlim([7 42]),ylim([9 53])
xlabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\DeltaT_{lat} (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
cc = corr(X(:,2),Y(:,2));
text(34,50,sprintf('r = %.2f',cc),'FontName','Arial','FontSize',13,'FontWeight','bold','VerticalAlignment','middle')
text(7.75,51,'A','FontName','Arial','FontSize',15,'FontWeight','bold','VerticalAlignment','middle')

% Panel B - Posterior gradients using posterior divisions
% subdivide using posterior divisions
ax2 = nexttile;
hold on, box on
for ii = 1:numel(fn)
    % 1 sigma bound
    y=[];
    for jj = 1:numel(CSpoidx.(fn{ii}))
        y = [y, LTG{CSpoidx.(fn{ii})(jj)}];
    end
    f = fill([Lat;flipud(Lat)],[prctile(y,16,2);...
        flipud(prctile(y,84,2))],cm(ii,:),...
        'FaceAlpha',.25,'EdgeColor','none');
    % Median line
    if fn(ii) == "tr"
        f.FaceColor = [.5, .5, .5];
        plot(Lat,median(y,2),'-','color',[.5, .5, .5],'LineWidth',3)
    else
        plot(Lat,median(y,2),'-','color',cm(ii,:),'LineWidth',3)
    end
end
% Add Text
x = 0;
y = [-15:-5:-40];
text(x,y(1),"PhanDA quantiles",'FontName','Arial','FontSize',13,...
    'HorizontalAlignment','center','Color','k','FontWeight','bold')
text(x,y(2),"Quantile #5 (hothouse)",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',cm(5,:))
text(x,y(3),"Quantile #4 (warmhouse)",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',cm(4,:))
text(x,y(4),"Quantile #3 (transitional)",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',[.5 .5 .5])
text(x,y(5),"Quantile #2 (coolhouse)",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',cm(2,:))
text(x,y(6),"Quantile #1 (coldhouse)",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',cm(1,:))
% Tidy plot
xlim([-90,90]);
ax2.XTick = [-90:30:90];
ax2.FontName = 'Arial';
ax2.FontSize = 11;
ylim([-50 46])
xlabel('Latitude','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['Surface air temperautre (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
text(-86,42,'B','FontName','Arial','FontSize',15,'FontWeight','bold','color','k')


% Panel C - Prior gradients using prior divisions
% subdivide using prior divisions
ax3 = nexttile;
hold on, box on
CSpridx.ih = find(GMSTprior<=Pprior(1));
CSpridx.ch = find(GMSTprior>Pprior(1) & GMSTprior<=Pprior(2));
CSpridx.tr = find(GMSTprior>Pprior(2) & GMSTprior<=Pprior(3));
CSpridx.gh = find(GMSTprior>Pprior(3) & GMSTprior<=Pprior(4));
CSpridx.hh = find(GMSTprior>Pprior(4));
for ii = 1:numel(fn)
    % 1 sigma bound
    y=[];
    for jj = 1:numel(CSpridx.(fn{ii}))
        y = [y, LTGprior(:,CSpridx.(fn{ii})(jj))];
    end
    if ~isempty(y)
        f = fill([Lat;flipud(Lat)],[prctile(y,16,2);...
            flipud(prctile(y,84,2))],cm(ii,:),...
            'FaceAlpha',.25,'EdgeColor','none');
        % Median line
        if fn{ii} == "tr"
            f.FaceColor = [.5, .5, .5];
            plot(Lat,median(y,2),'-','color',[.5, .5, .5],'LineWidth',3)
        else
            plot(Lat,median(y,2),'-','color',cm(ii,:),'LineWidth',3)
        end
    end
end
% Add Text
x = 0;
y = [-15:-5:-40];
text(x,y(1),"Model prior quantiles",'FontName','Arial','FontSize',13,...
    'HorizontalAlignment','center','Color','k','FontWeight','bold')
text(x,y(2),"Quantile #5",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',cm(5,:))
text(x,y(3),"Quantile #4",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',cm(4,:))
text(x,y(4),"Quantile #3",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',[.5 .5 .5])
text(x,y(5),"Quantile #2",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',cm(2,:))
text(x,y(6),"Quantile #1",'FontName','Arial','FontSize',11,...
    'HorizontalAlignment','center','Color',cm(1,:))
% Tidy plot
xlim([-90,90]);
ax3.XTick = [-90:30:90];
ax3.FontName = 'Arial';
ax3.FontSize = 11;
ylim([-50 46])
xlabel('Latitude','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['Surface air temperautre (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
text(-86,42,'C','FontName','Arial','FontSize',15,'FontWeight','bold','color','k')

% (e) Save figure
if savefig
    export_fig(gcf,[figdir,'/Supplemental/','SupFig_LTG.png'],'-p0.01','-m5')
end

