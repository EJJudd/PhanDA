%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Literature Comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure details
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
savefig = true;

% Select colormap
cm = hex2rgb({'#004F60';'#0a9396';'#ffb703';'#ca6702';'#9b2226'},1);

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

% Scotese 2021 (updated by van der Meer et al., 2022)
load("GMST_vanderMeer.mat")
% Valdes et al., 2021
load("GMST_Valdes.mat")
% Mills et al., 2021
load("GMST_Mills.mat")


%% MAKE FIGURE
figname = 'SupFig_LitCompPhanerozoic.png';
fig = figure('Position',[-1248,491,479*2,420],'Color','w'); 
hold on, box on
% Plot PhanDA GMST
GMSTmedian = cell2mat(cellfun(@(x) prctile(x,[5,16,50,84,95]), GMST, 'UniformOutput', false));
fill([GTS.Average;flipud(GTS.Average)],[GMSTmedian(:,1);...
    flipud(GMSTmedian(:,end))],'k','FaceColor',[.85 .85 .85],'EdgeColor','none');
fill([GTS.Average;flipud(GTS.Average)],...
    [GMSTmedian(:,2);flipud(GMSTmedian(:,end-1))],'k','FaceColor',...
    [.65 .65 .65],'EdgeColor','none');
p1 = plot(GTS.Average,GMSTmedian(:,3),'k-','LineWidth',2);
% Plot Scotese/van der Meer GMST
fill([GMST_vanderMeer.Age;flipud(GMST_vanderMeer.Age)],...
    [GMST_vanderMeer.GMST(:,1);flipud(GMST_vanderMeer.GMST(:,3))],...
    'k','FaceColor',cm(2,:),'EdgeColor','none','FaceAlpha',.25);
p2 = plot(GMST_vanderMeer.Age,GMST_vanderMeer.GMST(:,2),'-','LineWidth',2,'Color',cm(2,:));
% Plot Mills GMST
fill([GMST_Mills.Age;flipud(GMST_Mills.Age)],...
    [GMST_Mills.GMST(:,1);flipud(GMST_Mills.GMST(:,3))],...
    'k','FaceColor',cm(3,:),'EdgeColor','none','FaceAlpha',.25);
p3 = plot(GMST_Mills.Age,GMST_Mills.GMST(:,2),'-','LineWidth',2,'Color',cm(3,:));
p4 = plot(GMST_Valdes.Age,GMST_Valdes.Foster,'-','LineWidth',2,'Color',cm(4,:));
p5 = plot(GMST_Valdes.Age,GMST_Valdes.Smooth,'-','LineWidth',2,'Color',cm(5,:));
geologictimescale(0,GTS.LowerBoundary(end),...
    'normal','reverse',gca,'standard','stages','off',9,2)
leg = legend([p1, p2, p3, p4, p5], {'PhanDA', ...
    ['Scotese et al., 2021',newline,'(updated by van der Meer et al., 2022)'],...
    'Mills et al., 2021',...
    ['Valdes et al., 2021',newline,'(smoothed CO_2)'],...
    ['Valdes et al., 2021',newline,'(Foster et al., 2017 CO_2)']}, ...
    'Orientation','horizontal','Location','northoutside','FontName','Arial','FontSize',11);
set(gca,'FontName','Arial','FontSize',11)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['GMST (',char(176),'C)'],'FontName','Arial','FontSize',13,'FontWeight','bold')

% Save figure
if savefig
    export_fig(gcf,[figdir,'/',figname],'-p0.01','-m5')
end
