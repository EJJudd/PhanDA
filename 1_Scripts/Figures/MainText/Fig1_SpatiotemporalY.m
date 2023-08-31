%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Temporal and spatial distribution of Y %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
figname = 'Fig1_SpatiotemporalYvals.png';
savefig = true;

% Select colormap
PlotSpecs.Color = hex2rgb({'#004F60';'#0a9396';'#ffb703';'#ca6702';'#9b2226'},1);

% PART 1: LOAD DATA
% Directory details
assdate = '27Jul2023';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/InputWorkspaces/','YYe.mat'])
load([assdir,'/InputWorkspaces/','Data.mat'])
load([assdir,'/OutputWorkspaces/','Output.mat'],'TASprior','Ndata')
load("GTS2020_PETM.mat","GTS")
% Calculate prior GMST
GMSTprior = cell2mat(cellfun(@(x) latweightgmst(x)', TASprior, 'UniformOutput', false));

% PART 2 PRE-TREAT DATA
DArange = 1:size(GMSTprior);
%  (a) Assign bin range and means
    binrng = [GTS.UpperBoundary(DArange(1));GTS.LowerBoundary(DArange(1:end))];
    binmn = GTS.Average(DArange);
%  (b) Create list of unique proxy types
    fn = unique(Data.ProxyType);
    for ii = 1:numel(fn)
        Ny.(fn{ii}) = zeros(length(binrng)-1,1);
        Coords.(fn{ii}) = [];
    end
    sfn = fieldnames(UPD);
%  (c) For each time bin and proxy type, find the number of entries and
%      unique sampling sites
    for jj = 1:numel(sfn)
        fn = fieldnames(Y.(sfn{jj}));
    for ii = 1:numel(fn)
        % total entries:
        Ny.(fn{ii})(jj)=numel(Y.(sfn{jj}).(fn{ii}));
        Coords.(fn{ii}) = [Coords.(fn{ii}); [...
            UPD.(sfn{jj}).(fn{ii}).ModLat,UPD.(sfn{jj}).(fn{ii}).ModLon]];
    end
    end
%  (d) Combine arargonite and calcite isotope data into same proxy field
    Ny.d18ac = Ny.d18a+Ny.d18cforam+Ny.d18cmacro;
    Ny = rmfield(Ny, {'d18cforam','d18cmacro','d18a'});
    Coords.d18ac = [Coords.d18a;Coords.d18cforam;Coords.d18cmacro];
    Coords = rmfield(Coords, {'d18cforam','d18cmacro','d18a'});
    % List of fieldnames for plotting temporally
    tfn = fieldnames(Ny);
    tfn = [tfn(end);tfn(1:end-1)];
    %List of fieldnames for plotting spatially
    sfn = flipud(tfn);

%% PART 3: MAKE FIGURE
% (a) Define properties for each proxy and initialize figure
PlotSpecs.Proxy = fliplr(["U^{K'}_{37}","\fontsize{12}TEX\fontsize{8}86",...
    "Mg/Ca","\fontsize{12}δ^{18}O\fontsize{8}phosphate",...
    "\fontsize{12}δ^{18}O\fontsize{8}carbonate"]);
PlotSpecs.MarkerSize = [13:-2:5];
fig = figure('Units','inches','Position',[4 2 6.5*1.5 5.5*1.5],'color','w');
% (b) Plot assimilated N
ax1 = subplot('Position',[.07, .455, .76, .54]);
posY = zeros(numel(binmn),numel(tfn));
for ii = 1:numel(tfn)
    % for each stage, plot rectangle
    for jj = 13:numel(binrng)-1
        w = binrng(jj+1)-binrng(jj);
        h = Ny.(tfn{ii})(jj);
        posX = binrng(jj);
        posY(jj,ii+1)= posY(jj,ii) + h;
        rectangle('Position', [posX, posY(jj,ii), w, h], ...
            'EdgeColor', 'none', 'FaceColor', PlotSpecs.Color(ii,:));
    end
end
ax1.FontSize=11;ax1.FontName='Arial';
ylim([0 50])
ax1.YTick = [0:5:40];
geologictimescale(GTS.UpperBoundary(14),GTS.LowerBoundary(size(GMSTprior,1)),'normal','reverse',ax1,'standard','stages','off',10,1)    
ylabel('Assimilated values (N)','FontName','Arial','FontWeight','bold','FontSize',13)
Xtext = 470;
Ytext = linspace(44,30,5);
for ii = 1:numel(Ytext)
    text(ax1, Xtext,Ytext(ii),PlotSpecs.Proxy(ii),'FontName','Arial',...
        'FontSize',13,'FontWeight','bold','Color',PlotSpecs.Color(ii,:))
end
% (c) Map inset
axes('Position',[.24,.745,.375,.25]);
ax2 = worldmap('World');setm(ax2,'meridianlabel','off','parallellabel','off');gridm off
geoshow(shaperead('landareas', 'UseGeoCoords', true),'EdgeColor','none','FaceColor',[.75 .75 .75]);
for ii = 1:numel(sfn)
    Coords.(sfn{ii}) = unique(Coords.(sfn{ii}),'rows');
    plotm(Coords.(sfn{ii})(:,1),Coords.(sfn{ii})(:,2),'o','MarkerEdgeColor',...
        'none', 'MarkerSize',PlotSpecs.MarkerSize(ii),...
        'MarkerFaceColor',PlotSpecs.Color(end-ii+1,:))   
end
% (d) Plot assimilated N Zoom in
ax3 = subplot('Position', [.835, .455, .087, .54]);
posY = zeros(numel(binmn),numel(tfn));
yyaxis right
for ii = 1:numel(tfn)
    % for each stage, plot rectangle
    for jj = 1:13
        w = binrng(jj+1)-binrng(jj);
        h = Ny.(tfn{ii})(jj);
        posX = binrng(jj);
        posY(jj,ii+1)= posY(jj,ii) + h;
        rectangle('Position', [posX, posY(jj,ii), w, h], ...
            'EdgeColor', 'none', 'FaceColor', PlotSpecs.Color(ii,:));
    end
end
ax3.FontSize=11;ax3.FontName='Arial';
yt = ax3.YTick; yl = ax3.YLim;
yyaxis left
ax3.YTick=yt;ax3.YLim = yl;ax3.YTickLabel='';
geologictimescale(GTS.UpperBoundary(1),GTS.LowerBoundary(13),'normal','reverse',ax3,'standard','stages','off',10,1)    
yyaxis right
ax3.YTick = [0:25:150]; ax3.YTickLabel = [0:25:150];
ylim(ax3,[-15 152]); ax3.TickLength = [0.02,0.05];
ylabel('Assimilated values (N)','FontName','Arial','FontWeight','bold','FontSize',13)
% (e) Plot prior GMST
ax4 = axes('Position',[.07, .063, .76, .34]); hold on
f = fill([GTS.Average(1:91);flipud(GTS.Average(1:91))],...
    [min(GMSTprior,[],2);flipud(max(GMSTprior,[],2))],'k');
    f.FaceAlpha = .25; f.EdgeColor = 'none';
plot(GTS.Average(1:91),median(GMSTprior,2),'k-','LineWidth',2)
yl = ylim; yt = ax4.YTick;
geologictimescale(GTS.UpperBoundary(14),GTS.LowerBoundary(91),'normal',...
    'reverse',ax4,'standard','stages','off',5.5,1)    
ax4.FontSize=11;ax4.FontName='Arial'; box on
ylabel('Prior GMST (^oC)','FontName','Arial','FontWeight','bold','FontSize',13)
xlabel('Age (Ma)','FontName','Arial','FontWeight','bold','FontSize',13)
% (f) Prior zoom
ax5 = axes('Position',[.835, .063, .087, .34]); hold on
ylim(yl)
yyaxis right
f = fill([GTS.Average(1:91);flipud(GTS.Average(1:91))],...
    [min(GMSTprior,[],2);flipud(max(GMSTprior,[],2))],'k');
    f.FaceAlpha = .25; f.EdgeColor = 'none';
plot(GTS.Average(1:91),median(GMSTprior,2),'k-','LineWidth',2)
geologictimescale(GTS.UpperBoundary(1),GTS.LowerBoundary(13),'normal',...
    'reverse',ax5,'standard','stages','off',5.5,2)    
ax5.FontSize=11;ax5.FontName='Arial'; box on
yyaxis left, ax5.YTickLabel = [];ylim(ylim(ax4))
yyaxis right
ax5.YTick = yt; ax5.YTickLabel = yt;ax5.TickLength = [0.03,0.5];
ylabel('Prior GMST (^oC)','FontName','Arial','FontWeight','bold','FontSize',13)
% (g) Add labels
text(ax1,Xtext+10,47.5,'A','FontName','Arial',...
    'FontSize',15,'FontWeight','bold','Color','k')
text(ax1,375,47.5,'B','FontName','Arial',...
    'FontSize',15,'FontWeight','bold','Color','k')
text(ax4,Xtext+10,28,'C','FontName','Arial',...
    'FontSize',15,'FontWeight','bold','Color','k')
% (h) Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end
