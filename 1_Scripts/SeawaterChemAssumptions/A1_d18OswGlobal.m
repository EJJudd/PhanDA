%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  CONSTRUCTING THE TEMPORAL  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  OXYGEN ISOTOPE SEAWATER CURVE  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 01/21 (E. Judd)
% Last updated: 06/23 (E. Judd)

% Script to construct a time series of the global d18Oseawater
% This script uses data from:
% (1) Westerhold et al. (2020) 
%     (https://doi.org/10.1126/science.aba6853) 
% (2) Defliese (2020)
%     (https://doi.org/10.1016/j.epsl.2020.116661)

% The first half of this script (PART 1 - 4) construct the seawater curve,
% and the second half (PANEL 1-6) produces the supplemental figure.

% DEFINE DIRECTORIES & READ IN DATA
datadir = '/Users/emilyjudd/Documents/PhanDA/4_NonGlobalFiles/SeawaterChemAssumptions';
savedir = '/Users/emilyjudd/Documents/PhanDA/3_GlobalFiles/PSMs/seawaterchem';
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
% Benthic foram data
d18b = readtable( [datadir,'/d18Obenthic_Westerhold.csv'] );
CenoTime = d18b.Time;
% Snoball Earth correction
Snowball.d18o = readtable([datadir,'/SnowballCorrection_Defliese.csv']);
Snowball.d18o = table2array(Snowball.d18o);
Snowball.age = Snowball.d18o(1,:)';
Snowball.d18o(1,:) = []; 
Snowball.d18o = permute(Snowball.d18o,[2 1]);
% Geologic time scale
load("GTS2020_PETM.mat","GTS");

% PART 1: SCALE TO BENTHIC STACK TO RECONSTRUCT LATE CENOZOIC GLOBAL
%         SEAWATER COMPOSITION
% Detrend data
CenoSw = d18b.ISOBENd18oLOESSsmooth-d18b.ISOBENd18oLOESSsmooth(1);
% Scale data such that the value at the LGM is 1 permil 
% (Schrag et al., 1996)
CenoSw = CenoSw/(CenoSw(9));
% Fix any values below -1 permil at -1 permil to reflect ice-free 
% conditions (Shackleton and Kennett, 1975)
CenoSw(CenoSw<-1) = -1;

% PART 2: DIVIDE THE CENOZOIC RECORD BY PHASE
% Looking at the scaled record, it appears as if there are 3 distinct
% phases on glaciation: the "onset" (Rupelian - Langhian), the "body" 
% (Serravallian - Piacenzian) and the "peak" (Gelasian - Holocene)
Onsetidx = CenoTime > GTS.UpperBoundary(GTS.Stage == "Langhian") & ...
    CenoTime < GTS.LowerBoundary(GTS.Stage == "Rupelian");
Bodyidx = CenoTime > GTS.UpperBoundary(GTS.Stage == "Piacenzian") & ...
    CenoTime < GTS.LowerBoundary(GTS.Stage == "Serravallian");
Peakidx = CenoTime < GTS.LowerBoundary(GTS.Stage == "Gelasian");

% PART 3: ASSIGN SEAWATER VALUES FOR THE LATE CENOZOIC, EARLY PALEOZOIC AND
%         LATE PALEOZOIC ICEHOUSES BASED ON THE THREE PHASES
% We'll assign d18O values at a 1 myr time step, and then after we'll bin
% them by stage
% First, let's resample to get a distribution of values for each phase
N=200;
Peak(1,1:N) = datasample(CenoSw(Peakidx),N);
Body(1,1:N) = datasample(CenoSw(Bodyidx),N);
Onset(1,1:N) = datasample(CenoSw(Onsetidx),N);
PhanTime = [.5:ceil(max(GTS.LowerBoundary))]';
PhanSwMy = ones(numel(PhanTime),N)*-1;
Flag = zeros(numel(PhanTime),1);

% PART 3A: EARLY PALEOZOIC
% Based on: Finnegan et al. (2011), Boucot et al. (2013), and references
% therein
%Onset
PhanSwMy(PhanTime<=448.5&PhanTime>446.5,:) = [Onset;Onset];
Flag(PhanTime<=448.5&PhanTime>446.5,:) = [1;1];
%UpBody
PhanSwMy(PhanTime<=446.5&PhanTime>444.5,:) = [Body; Body];
Flag(PhanTime<=446.5&PhanTime>444.5,:) = [2; 2];
%Peak
PhanSwMy(PhanTime==444.5,:) = Peak;
Flag(PhanTime==444.5,:) = 3;
%DownBody
PhanSwMy(PhanTime==443.5,:) = Body;
Flag(PhanTime==443.5,:) = 2;
%Offset
PhanSwMy(PhanTime==442.5,:) = Onset;
Flag(PhanTime==442.5,:) = 1;

% PART 3B: LATE PALEOZOIC
% Based on Montenez et al. (2022) and references therein.
%Onset (end Famennian - mid Serpukhovian)
PhanSwMy(PhanTime<=361.5&PhanTime>325.5,:) = repmat(Onset,36,1);
Flag(PhanTime<=361.5&PhanTime>325.5,:) = ones(36,1);
%UpBody (mid Serpukhovian - mid Gzhelian)
PhanSwMy(PhanTime<=325.5&PhanTime>302.5,:) = repmat(Body,23,1);
Flag(PhanTime<=325.5&PhanTime>302.5,:) = repmat(2,23,1);
%Peak (mid Gzhelian - mid Asselian)
PhanSwMy(PhanTime<=302.5&PhanTime>296.5,:) = repmat(Peak,6,1);
Flag(PhanTime<=302.5&PhanTime>296.5,:) = repmat(3,6,1);
%DownBody (mid Asselian - end Sakmarian)
PhanSwMy(PhanTime<=296.5&PhanTime>290.5,:) = repmat(Body,6,1);
Flag(PhanTime<=296.5&PhanTime>290.5,:) = repmat(2,6,1);
%Offset (end Sakmarian - mid Kungurian)
PhanSwMy(PhanTime<=290.5&PhanTime>278.5,:) = repmat(Onset,12,1);
Flag(PhanTime<=290.5&PhanTime>278.5,:) = ones(12,1);

% Convert from 1 myr time step to stage level
PhanGlobalSW.Age = GTS.Average;
PhanGlobalSW.GlobalSW = interp1(PhanTime,PhanSwMy,PhanGlobalSW.Age);

% PART 3C: LATE CENOZOIC
% For the late Cenozoic, resample from the avaialble values for each stage
Bins = discretize(CenoTime,GTS.UpperBoundary);
for ii = 1:24
    PhanGlobalSW.GlobalSW(ii,:) = datasample(CenoSw(Bins==ii),N);
end

% PART 4: ADD IN THE DEFLIESE (2020) SNOWBALL EARTH CORRECTION
% Convert to stage level
PhanGlobalSW.Snowball = interp1(Snowball.age,Snowball.d18o,PhanGlobalSW.Age);
% Adjust data such that the modern correction is 0 per mil
adj = mean(PhanGlobalSW.Snowball(1,:));
PhanGlobalSW.Snowball = PhanGlobalSW.Snowball-adj;

% PART 5: ADD IN THE VEIZER AND PROKOPH (2015) CORRECTION
PhanGlobalSW.Veizer = -3e-5 * GTS.Average.^2 + 4.6e-3 * GTS.Average;

save([savedir,'/Phanerozoicd18Ov6.mat'],"PhanGlobalSW")

%% SECOND HALF: Plot figure
fig = figure('Position',[750 1 645 720],'Color','w');
tiledlayout(4,6,'Padding','none','TileSpacing','none');
% PANEL 1: Timeseries of scaled benthic stack
nexttile([1,4]), hold on
ylim([-1.1,1.2]), yl = ylim;
cm = hex2rgb({'#0A9396','#B4BE65','#FFB703'},1);
rectangle('Position',[min(CenoTime(Onsetidx)),yl(1),...
    range(CenoTime(Onsetidx)),range(yl)],...
    'FaceColor',[cm(3,:),.45],'EdgeColor','none')
text(mean(CenoTime(Onsetidx)),1.1,'Onset','Color',cm(3,:),'FontWeight',...
    'bold','FontSize',11,'FontName','Arial','HorizontalAlignment','center')
rectangle('Position',[min(CenoTime(Bodyidx)),yl(1),...
    range(CenoTime(Bodyidx)),range(yl)],...
    'FaceColor',[cm(2,:),.45],'EdgeColor','none')
text(mean(CenoTime(Bodyidx)),1.1,'Body','Color',cm(2,:),'FontWeight',...
    'bold','FontSize',11,'FontName','Arial','HorizontalAlignment','center')
rectangle('Position',[min(CenoTime(Peakidx)),yl(1),...
    range(CenoTime(Peakidx)),range(yl)],...
    'FaceColor',[cm(1,:),.45],'EdgeColor','none')
text(mean(CenoTime(Peakidx)),1.1,'Peak','Color',cm(1,:),'FontWeight',...
    'bold','FontSize',11,'FontName','Arial','HorizontalAlignment','center')
text(36.5,1.025,'A','Color','k','FontWeight',...
    'bold','FontSize',15,'FontName','Arial','HorizontalAlignment','center')
plot(CenoTime,CenoSw,'k-')
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Priabonian"),...
    'normal','reverse',gca,'standard','stages','off',5,1)
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\fontsize{13}δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',12,'FontWeight','bold')
% PANEL 2: Histogram of phases
nexttile([1,2]), hold on, box on
histogram(CenoSw(Onsetidx),'FaceColor',cm(3,:),'BinWidth',.1,'FaceAlpha',.45)
histogram(CenoSw(Bodyidx),'FaceColor',cm(2,:),'BinWidth',.1,'FaceAlpha',.45)
histogram(CenoSw(Peakidx),'FaceColor',cm(1,:),'BinWidth',.1,'FaceAlpha',.45)
xlabel(['\fontsize{13}δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel('Frequency (N)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylim([0 3500])
text(-.975,3310,'B','Color','k','FontWeight',...
    'bold','FontSize',15,'FontName','Arial','HorizontalAlignment','center')
% PANEL 3: Early Paleozoic (by myr)
nexttile([1,3]), hold on, box on
idx = find(PhanTime>400 & Flag ~= 0);
w = 1;
for ii = 1:numel(idx)
    x = PhanTime(idx(ii))-.5;
    y = min(PhanSwMy(idx(ii),:));
    h = range(PhanSwMy(idx(ii),:));
    rectangle('Position',[x,y,w,h],'FaceColor',[cm(4-Flag(idx(ii)),:),.45],...
        'EdgeColor','none')
end
plot(PhanTime,median(PhanSwMy,2),'k-','LineWidth',2)
ylim([-1.1 1.1])
geologictimescale(GTS.UpperBoundary(GTS.Stage == "Rhuddanian"),...
    GTS.LowerBoundary(GTS.Stage=="Katian"),'normal','reverse',gca,...
    'standard','stages','off',5,1)
text(452.3,.915,'C','Color','k','FontWeight',...
    'bold','FontSize',15,'FontName','Arial','HorizontalAlignment','center')
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\fontsize{13}δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',12,'FontWeight','bold')
% PANEL 4: Late Paleozoic (by myr)
nexttile([1,3]), hold on, box on
idx = find(PhanTime<400 & PhanTime>200 & Flag ~= 0);
w = 1;
for ii = 1:numel(idx)
    x = PhanTime(idx(ii))-.5;
    y = min(PhanSwMy(idx(ii),:));
    h = range(PhanSwMy(idx(ii),:));
    rectangle('Position',[x,y,w,h],'FaceColor',[cm(4-Flag(idx(ii)),:),.45],...
        'EdgeColor','none')
end
plot(PhanTime,median(PhanSwMy,2),'k-','LineWidth',2)
ylim([-1.1 1.1])
geologictimescale(GTS.UpperBoundary(GTS.Stage == "Kungurian"),...
    GTS.LowerBoundary(GTS.Stage=="Famennian"),'normal','reverse',gca,...
    'standard','stages','off',5,1)
text(367.5,.915,'D','Color','k','FontWeight',...
    'bold','FontSize',15,'FontName','Arial','HorizontalAlignment','center')
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\fontsize{13}δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
% PANEL 5: Phanerozoic seawater (by stage)
nexttile([1,6]), hold on, box on
fill([PhanGlobalSW.Age;flipud(PhanGlobalSW.Age)],...
    [prctile(PhanGlobalSW.GlobalSW,5,2);flipud(prctile(...
    PhanGlobalSW.GlobalSW,95,2))],'k','FaceAlpha',.25,'EdgeColor','none')
fill([PhanGlobalSW.Age;flipud(PhanGlobalSW.Age)],...
    [prctile(PhanGlobalSW.GlobalSW,16,2);flipud(prctile(...
    PhanGlobalSW.GlobalSW,84,2))],'k','FaceAlpha',.25,'EdgeColor','none')
plot(PhanGlobalSW.Age,median(PhanGlobalSW.GlobalSW,2),'color','k','LineWidth',2)
ylim([-1.1 1.1])
geologictimescale(0,GTS.LowerBoundary(91),'normal','reverse',gca,...
    'standard','stages','off',5,1)
text(478,.945,'E','Color','k','FontWeight',...
    'bold','FontSize',15,'FontName','Arial','HorizontalAlignment','center')
ylabel(['\fontsize{13}δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
% PANEL 6: Phanerozoic seawater with snoball correction (by stage)
nexttile([1,6]), hold on, box on
fill([PhanGlobalSW.Age;flipud(PhanGlobalSW.Age)],...
    [prctile(PhanGlobalSW.GlobalSW+mean(PhanGlobalSW.Snowball,2),5,2);flipud(prctile(...
    PhanGlobalSW.GlobalSW+mean(PhanGlobalSW.Snowball,2),95,2))],'k','FaceAlpha',.25,'EdgeColor','none')
fill([PhanGlobalSW.Age;flipud(PhanGlobalSW.Age)],...
    [prctile(PhanGlobalSW.GlobalSW+mean(PhanGlobalSW.Snowball,2),16,2);flipud(prctile(...
    PhanGlobalSW.GlobalSW+mean(PhanGlobalSW.Snowball,2),84,2))],'k','FaceAlpha',.25,'EdgeColor','none')
plot(PhanGlobalSW.Age,median(PhanGlobalSW.GlobalSW,2)+mean(PhanGlobalSW.Snowball,2),'color','k','LineWidth',2)
ylim([-2.45 1.1])
geologictimescale(0,GTS.LowerBoundary(91),'normal','reverse',gca,...
    'standard','stages','off',5,1)
text(478,.825,'F','Color','k','FontWeight',...
    'bold','FontSize',15,'FontName','Arial','HorizontalAlignment','center')
xlabel('Age (Ma)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\fontsize{13}δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
plot(xlim,[-1,-1],'k--')
box on
export_fig(fig,[figdir,'/SupFig_d18Global.png'],'-p0.01','-m5')
















