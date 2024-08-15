%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  Screening Protocol  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures';
savefig = true;
figname = 'SupFig_Screening.png';

% Select colormap
cm = hex2rgb({'#004F60','#0A9396','#c6c6c6','#CA6702','#9B2226'},1);

% PART 1: LOAD DATA
% Directory details
iCloud = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs';
ensdir = [iCloud,'/ModelOutputs/AssimilationFiles_beta'];
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
load([assdir,'/InputWorkspaces/YYe.mat'])
load("GTS2020_PETM.mat","GTS")
load("ExpNo.mat", "ExpNo")
load("HadCM3Coordinates.mat","Lat")
% Load ens
enfilename = [ensdir,'/HadCM3_all.ens'];
Ens = ensemble(enfilename);
EnsMeta = Ens.metadata;
dims = string(EnsMeta.ensembleDimensions{1});
ExpNoAll = EnsMeta.ensemble{1,1}{1, dims == "expno"};
% Select stage for assimilation
stageidx = 49;
stagename = GTS.Stage(stageidx);
stagelab = sprintf("S%2.0f_%s",stageidx,stagename);
% Isolate stage ensemble
expidx = [ExpNo.ExpNo1(stageidx), ExpNo.ExpNo2(stageidx)];
ensidx = any(ExpNoAll == expidx,2);
ensuse = Ens.useMembers(ensidx);
ensusemeta = ensuse.metadata;
ensusemat = load(ensuse);
exprunall = ensusemeta.ensemble{1,1}{1, dims == "exprun"};
% Define R values    
Rvals.d18cforam = [1e-3,1e-2,1e-1];
Rvals.d18acmacro = [5e-2,.275,5e-1];
Rvals.d18p = [5e-2,.275,5e-1; ...
              .1,.5,.9];
Rvals.tex = [1e-4,1e-3,1e-2];
Rvals.uk = [2.5e-5,2.5e-4,2.5e-3];
Rvals.mg = [5e-3,5e-2,5e-1; ...
                .1 .5 .9];
% Define assimilation options
Rmethod = "medium";
swcorr = "on";
pHcorr = "rec";
covanalysis = false;
% Calculate prior
tasprior = kel2cel(ensusemeta.regrid("tas",ensusemat,'order',["lat";"lon"]));
% Assemble Y, Ye, R
[Yuse, Yeuse, Ruse, proxytype, paleolat, paleolon] = assembleYYeR( UPD.(stagelab), ...
        Y.(stagelab), Ye.(stagelab), Rvals, Assumptions.(stagelab), ...
        pHcorr, swcorr, Rmethod, stageidx);
[Index, exprun, exps, explabs] = indexensemble(ensusemeta,stageidx);
% Run assimilation
[gmst, ltg, ltgsst, taspost, ~, Index] = runkf(ensusemeta, ensusemat, ...
    Yuse, Ruse, Yeuse, Index, "tas", covanalysis);

%% Plot figure
cm = flipud(customcolormap(linspace(0,1,8),{'#59475C','#005f73','#0a9396','#b4be65','#ffb703','#ca6702','#bb3e03','#9b2226'},8));
fig = figure('Position',[-1872 1 977 679],'Color','w'); 
t = tiledlayout(7,7,'Padding','none','TileSpacing','none');

% (1) Plot the individual configurations
tilenum = [3:7,10:14,17:21,24:28,31:35,38:42,45:49];
for ii = 1:size(Index,1)
    postt = nexttile(tilenum(ii)); hold on, box on
    plotltg(exprun{ii},exps{ii},exps{1},ltg(:,Index{ii,2}),...
        gmst(Index{ii,2}),explabs{ii},Index{ii,3},cm,postt,false)
    postt.FontSize = 10;
end

% (2) Allocate remaining axes
titlet = nexttile([1,2]);
priort = nexttile([3,2]); hold on, box on
postt = nexttile([3,2]); hold on, box on

% (3) Plot the prior
for ii = 1:numel(exps{1})
    plot(priort,Lat,squeeze(mean(tasprior(:,:,exprunall==exps{1}(ii)),2)),'-',...
        'color',cm(ii,:))
end
yl = ylim(priort);
y = linspace(yl(1)+.1*range(yl),yl(2)-range(yl)/2,numel(exps{1})+4);
for ii = 1:numel(exps{1})    
    text(priort,0,y(ii),exps{1}(ii),'HorizontalAlignment','center','FontWeight','bold',...
        'Color',cm(ii,:),'FontSize',11)
end
plot(priort,Lat,squeeze(mean(median(tasprior,3),2)),'k--','LineWidth',4)
text(priort,0,y(end-3),'mean of all','HorizontalAlignment','center','FontWeight','bold',...
        'Color','k','FontSize',11)
text(priort,0,y(end),"Prior",'HorizontalAlignment','center','FontSize',13,'FontWeight','bold')
text(priort,0,y(end-1),sprintf('(GMST = %.0f%sC)',median(latweightgmst(tasprior)),char(176)),...
    'HorizontalAlignment','center','FontSize',13,'FontWeight','normal')
priort.FontSize = 10; xlim(priort,[-90 90]), priort.XTick = [-90:30:90];

% (4) Plot posterior ltg
l = [];
g=[];
for ii = 1:size(Index,1)
    if Index{ii,3} == 1
    for jj = 1:numel(exps{ii})
        color = cm(exps{1} == exps{ii}(jj),:);
        plot(postt,Lat,ltg(:,Index{ii,2}(exprun{ii}==exps{ii}(jj))),'-','color',color)
        l = [l,ltg(:,Index{ii,2}(exprun{ii}==exps{ii}(jj)))];
        g = [g;gmst(Index{ii,2}(exprun{ii}==exps{ii}(jj)))];
    end
    end
end
plot(postt,Lat,median(l,2),'k--','LineWidth',4)
postt.FontSize = 11;
xlim([-90 90])
postt.XTick = [-90:30:90];
text(postt,0,y(end),"Accepted posteriors",'HorizontalAlignment','center','FontSize',13,'FontWeight','bold')
text(postt,0,y(end-1),sprintf('(GMST = %.0f%sC)',median(g),char(176)),'HorizontalAlignment','center','FontSize',13,'FontWeight','normal')

% (5) Title
nlines = 7;
x = linspace(.9,0,nlines);
text(titlet,.5,x(1),"Screening protocol example",'FontSize',18,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle')
text(titlet,.5,x(3),sprintf("%s (%.1f-%.1f Ma)",stagename,...
    GTS.LowerBoundary(stageidx),GTS.UpperBoundary(stageidx)),'FontSize',18,...
    'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','middle')
text(titlet,.5,x(5),"Snowball correction: on",'FontSize',13,...
    'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','middle')
text(titlet,.5,x(6),['pH correction: CO','\fontsize{7}2', '\fontsize{13} proxy pH'],'FontSize',13,...
    'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','middle')
text(titlet,.5,x(7),"R value: medium",'FontSize',13,...
    'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','middle')
titlet.Visible = 'off';

xlabel(t,'Latitude','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(t,['Surface air temperature (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')

if savefig
    export_fig(gcf,[figdir,'/Supplemental/',figname],'-p0.01','-m5')
end



