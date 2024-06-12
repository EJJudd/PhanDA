%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   CONSTRUCTING THE LOCAL   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   OXYGEN ISOTOPE SEAWATER REGRESSION MODEL %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 01/21 (E. Judd)
% Last updated: 06/23 (E. Judd)

% Script to construct regression model to predict local d18Osw
% This script uses:
% (1) modern d18Osw data from Legrande and Schmidt (2006) 
%     (https://doi.org/10.1029/2006GL026011)
% (2) modern salinity from the WOA18 published by Zweng et al. (2019)
%     (https://www.ncei.noaa.gov/sites/default/files/2020-04/woa18_vol2.pdf)
% (3) preindustrial HadCM3L simulations from Tindall et al. (2010)
%     (https://doi.org/10.1016/j.epsl.2009.12.049)
% (4) Carboniferous iCESM runs from Macarewich et al. (2021)
%     (https://doi.org/10.1016/j.epsl.2021.116770)

% The first half of this script (PART 1 - 4) constructs the seawater
% regression model and groundtruths it in the modern and Carboniferous,
% and the second half (PANEL 1-6) produces the supplemental figure.


% PART 1: DEFINE DIRECTORIES & READ IN DATA
% Directories
datadir = '/Users/emilyjudd/Documents/PhanDA/4_NonGlobalFiles/SeawaterChemAssumptions';
savedir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Code/DataAssimilation/4_GlobalFiles/PSMs/seawaterchem';
figdir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
% Load modern climate vars
load([datadir,'/ModernClimateVars.mat'])
% Remove boreal high latitude data
lat_idx = lat>75;
sss_anomaly(lat_idx) = NaN;
precipln_anomaly(lat_idx) = NaN;
dist(lat_idx) = NaN;
sst_anomaly(lat_idx) = NaN;
d18w_obs(lat_idx) = NaN;
% Load Ocean basin tag
load([datadir,'/OceanTag.mat'])


% PART 2: RUN MODEL
%  Define constants and preallocate variables
    %Radius of Earth (in km)
    R = 6371; 
    %Grid size
    s = size(sss_anomaly);
    %Search tolerance (in km)
    searchtol = 1750;
    % Model array
    mdl = cell(s);
    % R2 matrix
    R2 = NaN(s);
    % Number of regressed cells
    SN = NaN(s);

%Iterater through regression
for ii = 1:s(1)
    tic
    for jj = 1:s(2)
        if ~isnan(sss_anomaly(ii,jj))
            if abs(Lat(ii,1))>55 || OceanTag(ii,jj) == "ep"
                searchtol = 1000;
            else
                searchtol = 2000;
            end
            D = reshape(distance(Lat(ii,jj),lon(ii,jj),Lat(:),lon(:),[R 0]),s);
            searchidx = find(D<=searchtol & ~isnan(sss_anomaly) & ...
                ~isnan(d18w_obs));
            while numel(searchidx) < 10
                searchtol = searchtol+250;
                searchidx = find(D<=searchtol & ~isnan(sss_anomaly) & ...
                    ~isnan(d18w_obs)); 
            end
            x = [sss_anomaly(searchidx),...
                abs(Lat(searchidx)),...
                dist(searchidx),...
                precipln_anomaly(searchidx),...
                sst_anomaly(searchidx)];
            y = d18w_obs(searchidx);
            mdl{ii,jj} = stepwiselm(x,y,'quadratic','Verbose',0);
            R2(ii,jj) = mdl{ii,jj}.Rsquared.Adjusted;  
            SN(ii,jj) = numel(searchidx);
        end
    end
    elapsed_time = toc;
    fprintf('%.0f/73 complete (%.0f seconds)\n',ii,elapsed_time)
end
save([savedir,'/SeawaterMdl.mat'],"precipln_anomaly","sss_anomaly","sst_anomaly","dist","d18w_obs","R2","mdl","lat","lon")

% PART 3: GROUNDTRUTH IN MODERN
clearvars -except *dir
printprogress = true;
AssumptionFiles.predictsw = load("SeawaterMdl.mat");
load([datadir,'/ModernClimateVars.mat'])
nonanidx = find(~isnan(sss_anomaly));
x.dist = dist(nonanidx);
x.lat = Lat(nonanidx);
x.precipln_anomaly = precipln_anomaly(nonanidx);
x.sss_anomaly = sss_anomaly(nonanidx);
x.sst_anomaly = sst_anomaly(nonanidx);
d18O_predict = predictd18Osw(x, AssumptionFiles, printprogress);
d18O_predict_mod = NaN(size(AssumptionFiles.predictsw.d18w_obs));
d18O_predict_mod(nonanidx) = d18O_predict;
save([datadir,'/ModernGroundtruth.mat'],'d18O_predict_mod')
fprintf('\n%.0f%% of predictions are within %s0.25 of the observed value\n',...
    sum(abs(d18O_predict_mod-d18w_obs)<0.25,'all')/...
    sum(~isnan(sst_anomaly),'all')*100, char(177))
fprintf('\nThe latitudinally weighted global mean offset is %.02f\n',...
    latweightgmst(abs(d18O_predict_mod-d18w_obs)))

% PART 4: GROUNDTRUTH IN CARBONIFEROUS
clearvars -except *dir
printprogress = true;
AssumptionFiles.predictsw = load("SeawaterMdl.mat");
load([datadir,'/CarbClimateVars.mat'])
nonanidx = find(~isnan(sss_anomaly));
x.lat = Lat(nonanidx);
x.sst_anomaly = sst_anomaly(nonanidx);
x.dist = dist(nonanidx);
x.sss_anomaly = sss_anomaly(nonanidx);
x.precipln_anomaly = precipln_anomaly(nonanidx);
% Run model and save output
d18O_predict = predictd18Osw(x, AssumptionFiles, printprogress);
d18O_predict_carb = NaN(size(Lat));
d18O_predict_carb(nonanidx) = d18O_predict;
save([datadir,'/CarboniferousGroundtruth.mat'],'d18O_predict_carb')
load([datadir,'/CarbCESM.mat'])
fprintf('\n%.0f%% of predictions are within %s0.25 of the observed value\n',...
    sum(abs(d18O_predict_carb-d18O_CESM_regridded)<0.25,'all')/...
    sum(~isnan(d18O_CESM_regridded+d18O_predict_carb),'all')*100, char(177))
fprintf('\n%.0f%% of predictions are within %s0.5 of the observed value\n',...
    sum(abs(d18O_predict_carb-d18O_CESM_regridded)<0.5,'all')/...
    sum(~isnan(d18O_CESM_regridded+d18O_predict_carb),'all')*100, char(177))
fprintf('\nThe latitudinally weighted global mean offset is %.02f\n',...
    latweightgmst(abs(d18O_predict_carb-d18O_CESM_regridded)))

% PART 5: GROUNDTRUTH IN EOCENE
clearvars -except *dir
printprogress = true;
AssumptionFiles.predictsw = load("SeawaterMdl.mat");
load([datadir,'/EoceneClimateVars.mat'])
nonanidx = find(~isnan(sss_anomaly));
x.lat = Lat(nonanidx);
x.sst_anomaly = sst_anomaly(nonanidx);
x.dist = dist(nonanidx);
x.sss_anomaly = sss_anomaly(nonanidx);
x.precipln_anomaly = precipln_anomaly(nonanidx);
% Run model and save output
d18O_predict = predictd18Osw(x, AssumptionFiles, printprogress);
d18O_predict_eocene = NaN(size(Lat));
d18O_predict_eocene(nonanidx) = d18O_predict;
% Subtract 1 (to account for ice-free conditions)
d18O_predict_eocene = d18O_predict_eocene - 1;
save([datadir,'/EoceneGroundtruth.mat'],'d18O_predict_eocene')
load([datadir,'/EoceneCESM.mat'])
fprintf('\n%.0f%% of predictions are within %s0.25 of the observed value\n',...
    sum(abs(d18O_predict_eocene-d18O_CESM_regridded)<0.25,'all')/...
    sum(~isnan(d18O_CESM_regridded+d18O_predict_eocene),'all')*100, char(177))
fprintf('\n%.0f%% of predictions are within %s0.5 of the observed value\n',...
    sum(abs(d18O_predict_eocene-d18O_CESM_regridded)<0.5,'all')/...
    sum(~isnan(d18O_CESM_regridded+d18O_predict_eocene),'all')*100, char(177))
fprintf('\n%.0f%% of sub-arctic predictions are within %s0.25 of the observed value\n',...
    sum(abs(d18O_predict_eocene(Lat(:,1)<66.5,:)-d18O_CESM_regridded(Lat(:,1)<66.5,:))<0.25,'all')/...
    sum(~isnan(d18O_CESM_regridded(Lat(:,1)<66.5,:)+d18O_predict_eocene(Lat(:,1)<66.5,:)),'all')*100, char(177))
fprintf('\n%.0f%% of sub-arctic predictions are within %s0.5 of the observed value\n',...
    sum(abs(d18O_predict_eocene(Lat(:,1)<66.5,:)-d18O_CESM_regridded(Lat(:,1)<66.5,:))<0.5,'all')/...
    sum(~isnan(d18O_CESM_regridded(Lat(:,1)<66.5,:)+d18O_predict_eocene(Lat(:,1)<66.5,:)),'all')*100, char(177))
fprintf('\nThe latitudinally weighted global mean offset is %.02f\n',...
    latweightgmst(abs(d18O_predict_eocene-d18O_CESM_regridded)))




%% SECOND HALF: FIGURE
clearvars -except *dir
load([datadir,'/ModernClimateVars.mat'])
load([datadir,'/ModernGroundtruth.mat'])
load("HadCM3Coordinates.mat")
[Lon,Lat] = meshgrid(Lon,Lat);
fig = figure('Position',[-1248,491,620,290],'Color','w');
tiledlayout(1,2,'Padding','none','TileSpacing','none');
cm = hex2rgb({'#018571','#2A9887','#52AA9D','#7BBDB3','#A4D0C9','#CCE2DF','#f5f5f5',...
    '#E8DCD1','#DBC4AC','#CEAB88','#C09263','#B37A3E','#A6611A'},1);
% PANEL 1: Observed salinity vs observed d18Osw
nexttile; hold on, box on
load([datadir,'/ModernSalinityAnomaly.mat'])
scatter(sssa_obs(:),d18w_obs(:),75,'filled','MarkerFaceColor',[.5 .5 .5],...
    'MarkerEdgeColor','k','MarkerFaceAlpha',.25,'MarkerEdgeAlpha',.1)
r = corr(sssa_obs(~isnan(sssa_obs)),d18w_obs(~isnan(d18w_obs)));
xlim([-35 25])
text(-33.5,1.5,'A','FontName','Arial','FontSize',15,'FontWeight','bold')
text(-33.5,1,sprintf('r = %.02f',r),'FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel('Observed SSS anomaly (PSU)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\fontsize{13}Observed δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
% PANEL 2: HacCM3 salinity vs observed d18Osw
nexttile; hold on, box on
scatter(sss_anomaly(:),d18w_obs(:),75,'filled','MarkerFaceColor',[.5 .5 .5],...
    'MarkerEdgeColor','k','MarkerFaceAlpha',.25,'MarkerEdgeAlpha',.1)
r = corr(sss_anomaly(~isnan(sss_anomaly)),d18w_obs(~isnan(d18w_obs)));
xlim([-35 25])
text(-33.5,1.5,'B','FontName','Arial','FontSize',15,'FontWeight','bold')
text(-33.5,1,sprintf('r = %.02f',r),'FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel('HadCM3L SSS anomaly (PSU)','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(['\fontsize{13}Observed δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontWeight','bold','Interpreter','tex')
export_fig(gcf,[figdir,'/SupFig_SSScomp.png'],'-p0.01','-m5')

% PANEL 3: Predicted - observed map
fig = figure('Position',[750 1 580 820],'Color','w');
tiledlayout(3,26,'Padding','none','TileSpacing','none');
nexttile([1,16]); 
ax = worldmap('World');mlabel off, plabel off, gridm off
pcolorm(Lat,Lon,d18O_predict_mod-d18w_obs)
caxis([-3.25 3.25])
colormap(cm)
shading interp, framem k
geoshow(shaperead('landareas', 'UseGeoCoords', true),'EdgeColor','k','FaceColor','none');
c = colorbar('south');
c.AxisLocation = 'out';
c.Ticks = [-3:1:3];
c.FontSize = 11;
xlabel(c,['\fontsize{13}Predicted - observed δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
% PANEL 4: Residuals
nexttile([1,10]), hold on, box on
histogram(d18O_predict_mod-d18w_obs,'BinWidth',.25,'FaceColor',[.5 .5 .5]);
ylabel('Frequency (N)','FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel({'\fontsize{13}Predicted - observed', ['\fontsize{13}δ^{18}O','\fontsize{7}sw','\fontsize{13}  (‰)']},'FontName','Arial','FontSize',13,'FontWeight','bold')
set(gca,'YTick',[0:250:1750])
text(-3.75,1720,'B','FontName','Arial','FontSize',15,'FontWeight','bold')
% PANEL 5: Predicted - CESM map
load([datadir,'/CarboniferousGroundtruth.mat'])
load([datadir,'/CarbCESM.mat'])
load([datadir,'/CarbClimateVars.mat'])
lsmask = sst_anomaly;lsmask(~isnan(lsmask)) = 0;lsmask(isnan(lsmask)) = 1;
nexttile([1,16]); 
ax = worldmap('World');mlabel off, plabel off, gridm off
pcolorm(Lat,Lon,d18O_predict_carb-d18O_CESM_regridded)
caxis([-3.25 3.25])
c.Ticks = [-3:1:3];
colormap(cm)
shading interp, framem k
contourm(Lat,Lon,lsmask,0,'k-','linewidth',1)
c = colorbar('south');
c.AxisLocation = 'out';
c.FontSize = 11;
xlabel(c,['\fontsize{13}Predicted - iCESM δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
% PANEL 6: Residuals
nexttile([1,10]), hold on, box on
histogram(d18O_predict_carb-d18O_CESM_regridded,'BinWidth',.25,'FaceColor',[.5 .5 .5]);
ylabel('Frequency (N)','FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel({'\fontsize{13}Predicted - iCESM', ['\fontsize{13}δ^{18}O','\fontsize{7}sw','\fontsize{13}  (‰)']},'FontName','Arial','FontSize',13,'FontWeight','bold')
text(-3.95,1910,'D','FontName','Arial','FontSize',15,'FontWeight','bold')

% PANEL 7: Predicted - CESM map
load([datadir,'/EoceneGroundtruth.mat'])
load([datadir,'/EoceneCESM.mat'])
load([datadir,'/EoceneClimateVars.mat'])
lsmask = sst_anomaly+d18O_CESM_regridded;lsmask(~isnan(lsmask)) = 0;lsmask(isnan(lsmask)) = 1;
nexttile([1,16]); 
ax = worldmap('World');mlabel off, plabel off, gridm off
pcolorm(Lat,Lon,d18O_predict_eocene-d18O_CESM_regridded)
caxis([-3.25 3.25])
c.Ticks = [-3:1:3];
colormap(cm)
shading interp, framem k
contourm(Lat,Lon,lsmask,0,'k-','linewidth',1)
c = colorbar('south');
c.AxisLocation = 'out';
c.FontSize = 11;
xlabel(c,['\fontsize{13}Predicted - iCESM δ^{18}O','\fontsize{7}sw', '\fontsize{13}  (‰)'],'FontName','Arial','FontSize',13,'FontWeight','bold')
% PANEL 8: Residuals
nexttile([1,10]), hold on, box on
histogram(d18O_predict_eocene-d18O_CESM_regridded,'BinWidth',.25,'FaceColor',[.5 .5 .5]);
ylabel('Frequency (N)','FontName','Arial','FontSize',13,'FontWeight','bold')
xlabel({'\fontsize{13}Predicted - iCESM', ['\fontsize{13}δ^{18}O','\fontsize{7}sw','\fontsize{13}  (‰)']},'FontName','Arial','FontSize',13,'FontWeight','bold')
text(-1.8,1150,'F','FontName','Arial','FontSize',15,'FontWeight','bold')


annotation('textbox', [.05 .975 .01 .01], 'String', 'A', 'FontName','Arial','FontSize',15,'FontWeight','bold','EdgeColor','none')
annotation('textbox', [.05 .64 .01 .01], 'String', 'C', 'FontName','Arial','FontSize',15,'FontWeight','bold','EdgeColor','none')
annotation('textbox', [.05 .3 .01 .01], 'String', 'E', 'FontName','Arial','FontSize',15,'FontWeight','bold','EdgeColor','none')


export_fig(fig,[figdir,'/SupFig_d18Local.png'],'-p0.01','-m5')





