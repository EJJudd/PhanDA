%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  CONSTRUCTING THE TEMPORAL  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CO2  CURVE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 01/21 (J. Tierney)
% Last updated: 06/23 (E. Judd)

% Script to construct a time series of the global pCO2 and seawater pH.
% This script uses data from:
% (1) Foster et al., (2017)
%     (https://doi.org/10.1038/ncomms14845)
% (2) Rae et al. (2021)
%     (https://doi.org/10.1146/annurev-earth-082420-063026)
% (3) Lenton et al. (2018)
%     (https://doi.org/10.1016/j.earscirev.2017.12.004)
% Additionally, this script has an option to estimate sw pH from the CO2 
% values using the CO2SYS carbonate system calculation function from 
% VanHeuven et al. (2011)
% (https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/CO2SYS/co2rprt.html)

% The first half of this script (PART 1 - 5) constructs the CO2 and
% seawater pH curves and the second half (PANEL 1-6) produces a
% supplemental figure.

% DEFINE DIRECTORIES & READ IN DATA
datadir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Code/DataAssimilation/5_NonGlobalFiles/SeawaterChemAssumptions';
savedir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Code/DataAssimilation/4_GlobalFiles/PSMs/seawaterchem';
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
fileNames = string({dir([datadir,'/co2*']).name});
load("GTS2020_PETM.mat")
types=["alkenone","boron","liverwort","stomata","phytane","paleosol","COPSE"]';
% PART 1: READ IN DATA
fileFull = strcat(datadir,'/',fileNames);
% Preallocate
age_all = []; age_low = []; age_high = []; 
co2_all = []; co2_low = []; co2_high = []; site_all = []; type_all = [];
% Our CO2 requirements are:
%   (a) CO2 estimates must be realistic (remove values < 150 ppm)
%   (b) CO2 estimates must have an upper error bound (this removes 4
%       stomatal values)
%   (c) Exclude any Cenozoic data not in Rae et al., 2021
%   (d) Exclude any suspect paleosol data (as indicated by the 'Quality
%       Flag' field
cutoff = 150;
ceno = GTS.LowerBoundary(GTS.Stage == "Danian");
paleo = GTS.LowerBoundary(GTS.Stage=="Induan");
for ii = 1:length(fileFull)
    datsNow = readtable(fileFull(ii));
    proxy = types(~cellfun(@isempty,regexp(fileFull(ii),types)));
    if proxy == "alkenone"
        idx = datsNow.DatabaseSource == "Rae" & datsNow.CO2 > cutoff;
    elseif proxy == "boron"
        idx = datsNow.DatabaseSource == "Rae" & datsNow.CO2 > cutoff;
    elseif contains(fileFull(ii),"paleosol")
        idx = datsNow.QualityFlag == 0 & datsNow.CO2 > cutoff & ...
            datsNow.Age > ceno & ~isnan(datsNow.CO2High1s);     
    else
        idx = datsNow.CO2 > cutoff & datsNow.Age > ceno & ~isnan(datsNow.CO2High1s);     
    end
    age_all = [age_all;datsNow.Age(idx)];
    age_low = [age_low;datsNow.AgeLower(idx)];
    age_high = [age_high;datsNow.AgeUpper(idx)];
    co2_all = [co2_all;datsNow.CO2(idx)];
    co2_low = [co2_low;datsNow.CO2Low1s(idx)];
    co2_high = [co2_high;datsNow.CO2High1s(idx)];
    site_all = [site_all;datsNow.Site(idx)];  
    type_all = [type_all;repmat(proxy,sum(idx),1)];
end
% Add an age uncertainty to Paleozoic data missing age bounds of 0.5 myrs
age_low(isnan(age_low)&age_all>paleo) = age_all(isnan(age_low)&age_all>paleo) + 0.5;
age_high(isnan(age_high)&age_all>paleo) = age_all(isnan(age_high)&age_all>paleo) - 0.5;
% Increase lower bound CO2 uncertainties that are unrealistic
co2_low(co2_low<100) = 100;

% PART 2: LOAD COPSE DATA
load([datadir,'/CopseCO2.mat']);
% Remove all but the earliest data
idx = CopseAge <= GTS.LowerBoundary(GTS.Stage=="Fortunian") & ...
    CopseAge >= GTS.UpperBoundary(GTS.Stage=="Darriwilian");
% COPSE uncertainty comes from different modeling scenarios and is fairly
% large in the Paleozoic. We cannot reproduce the plot in Lenton 2018
% exactly but to approximate, we assign a healthy 1-sigma uncertainty of 
% 700 ppm.
CopseUnc = 700;
% Pre-bin the COPSE data by stage
StageEdges = [0; GTS.LowerBoundary];
StageMiddle = GTS.Average;
[~,~,loc]=histcounts(CopseAge(idx),StageEdges);
CopseCO2 = accumarray(loc(:),CopseCO2(idx))./accumarray(loc(:),1);
age_all = [age_all;StageMiddle(unique(loc))];
age_high = [age_high;NaN(size(CopseCO2(~isnan(CopseCO2))))];
age_low = [age_low;NaN(size(CopseCO2(~isnan(CopseCO2))))];
co2_all = [co2_all; CopseCO2(~isnan(CopseCO2))];
co2_high= [co2_high; CopseCO2(~isnan(CopseCO2))+CopseUnc];
co2_low = [co2_low; CopseCO2(~isnan(CopseCO2))-CopseUnc];
type_all = [type_all; repmat("COPSE",size(CopseCO2(~isnan(CopseCO2))))];
site_all = [site_all; repmat("COPSE",size(CopseCO2(~isnan(CopseCO2))))];

% PART 3: RUN THE MONTE CARLO SIMULATION
Niters = 10000;
CO2_binned = NaN(length(StageMiddle),Niters);
idx = ~isnan(age_low);
for ii=1:Niters
    age_now = age_all;
    age_now(idx) = unifrnd(age_high(idx),age_low(idx));
    bins = discretize(age_now,StageEdges);
    ubin = unique(bins);
    for jj = 1:length(ubin)
        bin_now = bins == ubin(jj);
        % ID the sites and proxy types in the bin
        site_now = site_all(bin_now);
        type_now = type_all(bin_now);
        siteType = strcat(site_all(bin_now),type_all(bin_now));
        usites = unique(siteType,'rows');
        % Grab co2 and errors in bin, convert to log space
        co2_now = log(co2_all(bin_now));
        low_now = log(co2_all(bin_now))-log(co2_low(bin_now));
        high_now = log(co2_high(bin_now))-log(co2_all(bin_now));
        % Take the mean log error as a an approx of 1-sigma
        avg_err = mean([low_now,high_now],2);
        % Pre-average by type and site so as not to bias outcome
        sitemean = NaN(length(usites),1);
        sitestd = NaN(length(usites),1);
        for Nfinal = 1:length(usites)
            sidx = siteType == usites(Nfinal,:);
            sitemean(Nfinal) = sum(co2_now(sidx)./(avg_err(sidx).^2))./sum(1./(avg_err(sidx).^2));
            % taking the max error from each site (conservative)
            sitestd(Nfinal) = max(avg_err(sidx));
        end
        % Weighted mean for all
        mean_now = sum(sitemean./(sitestd.^2))./sum(1./(sitestd.^2));
        %error prop assuming independence
        std_now = sqrt(1./sum(1./(sitestd.^2)));
        % sample from the distribution
        CO2_binned(ubin(jj),ii) = exp(normrnd(mean_now,std_now));
    end
    if mod(ii,500) == 0
        fprintf('Completed %.0f of %.0f iterations\n',ii,Niters)
    end
end

% PART 4: REMOVE NANS AND RESIZE
Nfinal = 2500;
PhanerozoicCO2 = NaN(numel(StageMiddle),Nfinal);
for ii = 1:size(PhanerozoicCO2,1)
    if sum(isnan(CO2_binned(ii,:))) <= Nfinal
        PhanerozoicCO2(ii,:) = datasample(CO2_binned(ii,~isnan(CO2_binned(ii,:))),...
            Nfinal,'Replace',false);
    elseif sum(isnan(CO2_binned(ii,:))) > Nfinal && sum(isnan(CO2_binned(ii,:))) ~= Niters
        PhanerozoicCO2(ii,:) = datasample(CO2_binned(ii,~isnan(CO2_binned(ii,:))),...
            Nfinal,'Replace',true);    
    else
        PhanerozoicCO2(ii,1:Nfinal/2) = datasample(CO2_binned(ii+1,~isnan(CO2_binned(ii+1,:))),...
            Nfinal/2,'Replace',false);
        PhanerozoicCO2(ii,Nfinal/2+1:end) = datasample(CO2_binned(ii-1,~isnan(CO2_binned(ii-1,:))),...
            Nfinal/2,'Replace',false);
    end
end

save([savedir,'/PhanerozoicCO2v9.mat'],"PhanerozoicCO2")

% PART 5: CALCULATE PH USING CO2SYS
% vectorize binned CO2 values
co2Vec = PhanerozoicCO2(:);
% Use modern distribution of alk.
alk = normrnd(2300,100,length(co2Vec),1);
% Random draws for T and S
T = unifrnd(10,35,length(co2Vec),1);
% Salinity based on modern distribution
S = normrnd(34,2,length(co2Vec),1);
% Define constants for carbonate system calculations
SCALE  = 1; % Total pH scale
K1K2   = 4; % K1 K2 from Mehrbach (1973), refit by Dickson & Millero (1987)
KSO4   = 1; % KSO4 from Dickson (1990) 
KF     = 2; % KHF from Perez & Fraga (1987)
BOR    = 2; % Borate-to-salinity ratio from Lee et al (2010)
% CO2SYS call. 
% Parameter 1 is alkalinity (1). Parameter 2 is CO2 (4). S and T come next,
% then Tout (and all outher 'out' values) is NaN. Surface pressure is 0.
% Silicate, phosphate, ammonium, H2S just set to arbitrary low values. 
A = CO2SYS(alk,co2Vec,1,4,S,T,nan,0,nan,1,1,0,0,SCALE,K1K2,KSO4,KF,BOR);
PhanerozoicpH = A(:,3);
PhanerozoicpH = reshape(PhanerozoicpH,size(PhanerozoicCO2));

save([savedir,'/PhanerozoicpHv6.mat'],"PhanerozoicpH")

%% SECOND HALF: Plot figure
fig = figure('Position',[750 1 645 820],'Color','w');
tiledlayout(4,1,'Padding','none','TileSpacing','compact');
cm_grey = customcolormap(linspace(0,1,2),{'#6A6A6A','#FCFCFC'},75);
cm_grey = [cm_grey;0,0,0;0,0,0;flipud(cm_grey)];
% Color for extinctions
cm_col = hex2rgb(['#004F60';'#0a9396';'#b4be65';'#ffb703';'#ca6702';'#9B2226';'#000000'],1);

% PANEL 1: CO2 Proxy data
ax1 = nexttile; hold on
fn = unique(type_all);
fn = [fn(2:end);fn(1)];
fnorder = [2,1,4,5,6,3,7];
collab = [];
for ii = 1:numel(fn)
    collab = strcat(collab,'\color[rgb]{',num2str(cm_col(ii,:)),'}',fn(ii),{'   '});
end
for ii = 1:numel(fn)
    idx = type_all == fn(fnorder(ii));
    errorbar(age_all(idx),log(co2_all(idx)),log(co2_all(idx)) - log(co2_low(idx)),...
        log(co2_high(idx))-log(co2_all(idx)),age_all(idx)-age_low(idx),...
        age_high(idx)-age_all(idx),'o','Color',cm_col(fnorder(ii),:),'LineWidth',.25,...
        'MarkerFaceColor','none','MarkerEdgeColor','none','CapSize',0)
end
for ii = 1:numel(fn)
    idx = type_all == fn(fnorder(ii));
    scatter(age_all(idx),log(co2_all(idx)),50,'filled','MarkerFaceColor',...
        cm_col(fnorder(ii),:),'MarkerEdgeColor','none','MarkerFaceAlpha',1)
end
ylim([4.5 9.4])
text(mean([0,GTS.LowerBoundary(GTS.Stage=="Tremadocian")]),9.1,collab,...
    'FontName','Arial','FontSize',13,'HorizontalAlignment','center')
ax1.YTick = log([100,200,500,1000,2000,5000]);
ax1.YAxis.MinorTick = 'on';
ax1.YAxis.MinorTickValues = log([100:100:1000,2000:1000:5000]);
ax1.YTickLabel = ["100","200","500","1000","2000","5000"];
ax1.FontSize = 11; ax1.FontName = 'Arial';
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Tremadocian"),'normal',...
    'reverse',ax1,'standard','stages','off',8,1)
ax1.YAxis(2).MinorTick = 'on';
ax1.YAxis(2).MinorTickValues = log([100:100:1000,2000:1000:5000]);
ylabel(['\fontsize{13}CO','\fontsize{7}2', '\fontsize{13}  (ppmv)'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
text(ax1,483.5,9.1,'A','FontSize',15,'FontWeight','bold','FontName','Arial')

% PANEL 2: CO2 CURVE
ax2 = nexttile; hold on
P = [5:1:95];
CO2 = prctile(PhanerozoicCO2,P,2);
[P,age_grid] = meshgrid(P,GTS.Average);
contourf(age_grid, log(CO2), P, 152, 'LineColor', 'none')
plot(GTS.Average,log(median(PhanerozoicCO2,2)),'k-','LineWidth',2)
ylim(log([75,7500]))
ax2.YTick = log([100,200,500,1000,2000,5000]);
ax2.YAxis.MinorTick = 'on';
ax2.YAxis.MinorTickValues = log([100:100:1000,2000:1000:5000]);
ax2.YTickLabel = ["100","200","500","1000","2000","5000"];
ax2.FontSize = 11; ax1.FontName = 'Arial';
colormap(cm_grey)
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Tremadocian"),'normal',...
    'reverse',ax2,'standard','stages','off',8,1)
ax2.YAxis(2).MinorTick = 'on';
ax2.YAxis(2).MinorTickValues = log([100:100:1000,2000:1000:5000]);
ylabel(['\fontsize{13}CO','\fontsize{7}2', '\fontsize{13}  (ppmv)'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
box on
text(ax2,483.5,8.65,'B','FontSize',15,'FontWeight','bold','FontName','Arial')

% PANEL 3: Foster Comparison
ax3 = nexttile; hold on
datadir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Code/DataAssimilation/5_NonGlobalFiles';
load([datadir,'/FosterCO2.mat'])
F = discretize(FosterCO2(:,1),StageEdges);
FosterCO2(FosterCO2<0) = 1;
FosterCO2_Stage = NaN(size(GTS,1),1);
for ii = 1:size(GTS,1)
    try
    FosterCO2_Stage(ii) = mean(FosterCO2(F==ii,4));
    end
end
fill([GTS.Average;flipud(GTS.Average)],...
    log([prctile(PhanerozoicCO2,2.5,2);flipud(prctile(PhanerozoicCO2,97.5,2))]),...
    'k','FaceAlpha',.25,'EdgeColor','none')
fill([GTS.Average;flipud(GTS.Average)],...
    log([prctile(PhanerozoicCO2,16,2);flipud(prctile(PhanerozoicCO2,84,2))]),...
    'k','FaceAlpha',.25,'EdgeColor','none')
plot(GTS.Average,log(median(PhanerozoicCO2,2)),'k-','LineWidth',2)
fill([FosterCO2(:,1);flipud(FosterCO2(:,1))],...
    log([FosterCO2(:,2);flipud(FosterCO2(:,6))]),...
    cm_col(1,:),'FaceAlpha',.25,'EdgeColor','none')
fill([FosterCO2(:,1);flipud(FosterCO2(:,1))],...
    log([FosterCO2(:,3);flipud(FosterCO2(:,5))]),...
    cm_col(1,:),'FaceAlpha',.25,'EdgeColor','none')
plot(FosterCO2(:,1),log(FosterCO2(:,4)),'-','LineWidth',2,'Color',cm_col(1,:))
plot(GTS.Average,log(FosterCO2_Stage),'--','LineWidth',2,'Color',cm_col(2,:))
ylim(log([75,7500]))
ax3.YTick = log([100,200,500,1000,2000,5000]);
ax3.YAxis.MinorTick = 'on';
ax3.YAxis.MinorTickValues = log([100:100:1000,2000:1000:5000]);
ax3.YTickLabel = ["100","200","500","1000","2000","5000"];
ax3.FontSize = 11; ax1.FontName = 'Arial';
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Tremadocian"),'normal',...
    'reverse',ax3,'standard','stages','off',8,1)
ax3.YAxis(2).MinorTick = 'on';
ax3.YAxis(2).MinorTickValues = log([100:100:1000,2000:1000:5000]);
ylabel(['\fontsize{13}CO','\fontsize{7}2', '\fontsize{13}  (ppmv)'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
box on
text(ax3,483.5,8.65,'C','FontSize',15,'FontWeight','bold','FontName','Arial')
text(ax3,5,log(5700),'PhanDA','FontSize',13,'FontName','Arial',...
    'HorizontalAlignment','right')
text(ax3,5,log(3900),'Foster et al., 2017','FontSize',13,'FontName',...
    'Arial','Color',cm_col(1,:),'HorizontalAlignment','right')
text(ax3,5,log(2650),'Foster et al., 2017','FontSize',13,'FontName',...
    'Arial','Color',cm_col(2,:),'HorizontalAlignment','right')
text(ax3,5,log(1975),'(stage averaged)','FontSize',11,'FontName','Arial',...
    'Color',cm_col(2,:),'HorizontalAlignment','right')


% PANEL 4: Seawater pH curve 
ax4 = nexttile; hold on
P = [5:1:95];
pH = prctile(PhanerozoicpH,P,2);
[P,age_grid] = meshgrid(P,GTS.Average);
contourf(age_grid, pH, P, 152, 'LineColor', 'none')
plot(GTS.Average,median(PhanerozoicpH,2),'k-','LineWidth',2)
ax4.FontSize = 11; ax1.FontName = 'Arial';
colormap(cm_grey)
geologictimescale(0,GTS.LowerBoundary(GTS.Stage=="Tremadocian"),'normal',...
    'reverse',ax4,'standard','stages','off',8,1)
ylabel(['\fontsize{13}pH','\fontsize{7}sw'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
xlabel('\fontsize{13}Age (Ma)','FontName','Arial','FontWeight','bold','Interpreter','tex')
box on
text(ax4,483.5,8.45,'D','FontSize',15,'FontWeight','bold','FontName','Arial')

export_fig(fig,[figdir,'/SupFig_GlobalCO2.png'],'-p0.01','-m5')


