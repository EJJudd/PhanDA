%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PhanDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
% PART 1: LOAD DATA
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","LTG","ItName","LTGsst","TASprior")
load([assdir,'/InputWorkspaces/','YYe.mat'],"Preferences")
load("GTS2020_PETM.mat","GTS")
load("HadCM3Coordinates.mat")
load("PhanerozoicCO2v9.mat", "PhanerozoicCO2")
load("lsmask.mat","lsmask");
lsmask = struct2cell(lsmask);

% PART 2: PRE-TREAT DATA
% Select iterations to use
pHCorr = ["ens","rec"];
swCorr = ["snowball","off"];
rMeth = ["low","medium","high"];
idx = contains(ItName,strcat("phCorr = ",pHCorr)) & ...
    contains(ItName,strcat("SeawaterCorr = ",string(swCorr))) & ...
    contains(ItName,strcat("Rmethod = ",rMeth));
GMST = combineruns(GMST,idx,1);
LTG = combineruns(LTG,idx,2);
LTGsst = combineruns(LTGsst,idx,2);

% Revise LTG, GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
LTG = combinestages(LTG,"LTG",Preferences,endsize,2);
LTGsst = combinestages(LTGsst,"LTGsst",Preferences,endsize,2);
PhanerozoicCO2 = combinestages(PhanerozoicCO2,"CO2",Preferences,endsize,2);
TASprior = combinestages(TASprior,"TASprior",Preferences,endsize,3);
lsmask = combinestages(lsmask,"lsmask",Preferences,endsize,3);

% Define the eras
cenozoic = find(GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Selandian/Danian"));
mesozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Selandian/Danian") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Induan"));
paleozoic = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<=GTS.LowerBoundary(GTS.Stage=="Tremadocian"));
paleonocopse = find(GTS.LowerBoundary>GTS.LowerBoundary(GTS.Stage=="Induan") & ...
    GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));
allnocopse = find(GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian"));
allnocopsemeso = setdiff(find(GTS.LowerBoundary<GTS.LowerBoundary(GTS.Stage=="Darriwilian")),mesozoic);

% Define the regions
poles = abs(Lat)>=66.5;
tropics = abs(Lat)<=23.5;

% Calculate percentiles
p = [16,50,84];
pgmst = cell2mat(cellfun(@(x) prctile(x,p), GMST, 'UniformOutput', false));
pco2 = prctile(log2(PhanerozoicCO2),p,2);

% subdivide median GMST by climate state quantile
q = prctile(pgmst(:,p==50),[20:20:100]);
CSidx.ih = find(pgmst(:,p==50)<=q(1));
CSidx.ch = find(pgmst(:,p==50)>q(1) & pgmst(:,p==50)<=q(2));
CSidx.tr = find(pgmst(:,p==50)>q(2) & pgmst(:,p==50)<=q(3));
CSidx.gh = find(pgmst(:,p==50)>q(3) & pgmst(:,p==50)<=q(4));
CSidx.hh = find(pgmst(:,p==50)>q(4));

% Proportion of time spent in each state in each era
fn = fieldnames(CSidx);
for ii = 1:numel(fn)
	CStime.(fn{ii}) = sum(GTS.LowerBoundary(CSidx.(fn{ii})) - ...
        GTS.UpperBoundary(CSidx.(fn{ii})))./GTS.LowerBoundary(end)*100;
end

% Calculate ACS & Correlation coefficients
% All
[~,mA,~,smA] = york_fit(pco2(allnocopse,p==50)',pgmst(allnocopse,p==50)',...
    (range(pco2(allnocopse,:),2)./2)',(range(pgmst(allnocopse,:),2)./2)', 0);
[pvalA,rA,~,cv] = ebisuzaki(pco2(allnocopse,p==50),pgmst(allnocopse,p==50));
% the absolute value of r is greater than the critical value and the
% p-value is <<.05, so the correlation is significant
% Cenozoic
[~,mC,~,smC] = york_fit(pco2(cenozoic,p==50)',pgmst(cenozoic,p==50)',...
    (range(pco2(cenozoic,:),2)./2)',(range(pgmst(cenozoic,:),2)./2)', 0);
[pvalC,rC,~,cv] = ebisuzaki(pco2(cenozoic,p==50),pgmst(cenozoic,p==50));
% the absolute value of r is greater than the critical value and the
% p-value is <<.05, so the correlation is significant
% Mesozoic
[pvalM,rM,~,cv] = ebisuzaki(pco2(mesozoic,p==50),pgmst(mesozoic,p==50));
% the absolute value of r is greater than the critical value and the
% p-value is 0.175, so the correlation is not significant
% Paleozoic
[~,mP,~,smP] = york_fit(pco2(paleonocopse,p==50)',pgmst(paleonocopse,p==50)',...
    (range(pco2(paleonocopse,:),2)./2)',(range(pgmst(paleonocopse,:),2)./2)', 0);
[pvalP,rP,~,cv] = ebisuzaki(pco2(paleonocopse,p==50),pgmst(paleonocopse,p==50));
% the absolute value of r is greater than the critical value and the
% p-value is 0.175, so the correlation is not significant
% the absolute value of r is less than the critical value and the
% p-value is <<.05, so the correlation is significant
% Cenozoic and Paleozoic combined
[~,mCP,~,smCP] = york_fit(pco2([paleonocopse;cenozoic],p==50)',...
    pgmst([paleonocopse;cenozoic],p==50)',(range(pco2([paleonocopse;cenozoic],:),2)./2)',...
    (range(pgmst([paleonocopse;cenozoic],:),2)./2)', 0);


% Time weight GMST
[TWave, TWstd] = timeweight(pgmst(:,p==50),[GTS.UpperBoundary(1);GTS.LowerBoundary]);

% Time weight mean prior GMST
GMSTprior = cell2mat(cellfun(@(x) latweightgmst(x), TASprior, ...
    'UniformOutput', false));
[TWaveprior, TWstdprior] = timeweight( cell2mat( cellfun(@(x) ...
    median(latweightgmst(x)), TASprior,'UniformOutput', false)),...
    [GTS.UpperBoundary(1);GTS.LowerBoundary]);

% Cenozoic literature comparison
[~, ~, ~, ~, PhanDA, Lit] = loadcenolitdata(GMST,GTS);

% Benthic stack
d18b = readtable('/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Data/SupplementalData/d18Obenthic_Westerhold.csv');
d18bStage = NaN(size(cenozoic));
for ii = 1:numel(cenozoic)
    idx = d18b.Time>=GTS.UpperBoundary(ii) & ...
        d18b.Time<GTS.LowerBoundary(ii);
    d18bStage(ii) = median(d18b.ISOBENd18oLOESSsmooth(idx));
end
[pvalBS,rBS,~,cv] = ebisuzaki(pgmst(cenozoic,p==50),d18bStage);
% the absolute value of r is greater than the critical value and the
% p-value is <<.05, so the correlation is significant

% SAT and SST LTG
ltgsat = NaN(size(GMST));
ltgsst = NaN(size(GMST));
tropsat = NaN(size(GMST));
tropsst = NaN(size(GMST));
meandim = 1;
for ii = 1:numel(GMST)
    xt = LTG{ii}; xt(~tropics,:) = NaN;
    xp = LTG{ii}; xp(~poles,:) = NaN;
    ltgsat(ii) = median(latweightgmst(xt,meandim)-latweightgmst(xp,meandim));
    tropsat(ii) = median(latweightgmst(xt));
    xt = LTGsst{ii}; xt(~tropics,:) = NaN;
    xp = LTGsst{ii}; xp(~poles,:) = NaN;
    ltgsst(ii) = median(latweightgmst(xt,meandim)-latweightgmst(xp,meandim));
    tropsst(ii) = median(latweightgmst(xt));
end
[pvalLTG,rLTG,~,cv] = ebisuzaki(pgmst(:,p==50),ltgsat);
% the absolute value of r is greater than the critical value and the
% p-value is <<.05, so the correlation is significant

% Forcing calculations
albedo = 0.29;
[~, dFtot, dFco2, dFsun, dFsa] = ...
    calcforcings(GTS.Average,PhanerozoicCO2,albedo,lsmask);


%% Statistics (in the order they appear in the paper
% Abstract:
% (1) Min and max GMST
fprintf('GMST fluctuates between %.0f and %.0fC\n', min(pgmst(:,p==50)), ...
    max(pgmst(:,p==50)));
% (2) Apparent climate sensitivity
fprintf('The ACS is %.1f\n', mA);

% The PhanDA reconstruction:
% (1) time-weighted mean GMST across all assimilated stages
fprintf("\nThe time-weighted mean GMST is: %.0f %s %.0f (1 sigma)\n", ...
    TWave, char(177), TWstd)
% (2) mean + std of prior GMST
fprintf("The prior mean GMST is: %.0f %s %.0f (1 sigma)\n",...
    mean(GMSTprior), char(177), std(GMSTprior))
fprintf("The time-weighted prior GMST is: %.0f %s %.0f (1 sigma)\n", ...
    TWaveprior, char(177), TWstdprior)
% (3) Cenozoic literature comparison
fprintf("The correlation coefficient between PhanDA and the Cenozoic literature is: %.2f\n",...
    corr(PhanDA(:,2),Lit(:,2)))

% Phanerozoic temperatures:
% (1) min/max GMST and stage
fprintf("\nThe minimum reconstructed GMST of %.0fC occurred in the %s.\n",...
    min(pgmst(:,p==50)),GTS.Stage(pgmst(:,p==50) == min(pgmst(:,p==50))))
fprintf("The maximum reconstructed GMST of %.0fC occurred in the %s.\n",...
    max(pgmst(:,p==50)),GTS.Stage(pgmst(:,p==50) == max(pgmst(:,p==50))))
% (2) correlation coefficient between GMST and the benthic stack
fprintf("The correlation coefficient between PhanDA and the benthic stack is: %.2f\n",...
    rBS)


% Phanerozoic climate states and latitudinal gradients
% (1) percent of time spent in icehouse and greenhouse states
fprintf("%.0f%s of the Phanerozoic was spent in a greenhouse state\n",...
    CStime.gh+CStime.hh,"%")
fprintf("%.0f%s of the Phanerozoic was spent in an icehouse state\n",...
    CStime.ch+CStime.ih,"%")
fprintf("%.0f%s of the Phanerozoic was spent in a coldhouse state\n",...
    CStime.ih,"%")
% (2) correlation coefficient between GMST and LTG
fprintf("The correlation coefficient between GMST and the LTG is: %.2f (p-value = %.04f)\n",...
    rLTG,pvalLTG)
% (3) coldhouse LTG range
fprintf("The LTG ranges between %.0f and %.0fC during coldhouse intervals\n",...
    min(ltgsat(CSidx.ih)),max(ltgsat(CSidx.ih)))
% (4) hothouse LTG range
fprintf("The LTG ranges between %.0f and %.0fC during hothouse intervals\n",...
    min(ltgsat(CSidx.hh)),max(ltgsat(CSidx.hh)))
% (5) icehouse LTG range
fprintf("The LTG ranges between %.0f and %.0fC during icehouse intervals\n",...
    min(ltgsat([CSidx.ch;CSidx.ih])),max(ltgsat([CSidx.ch;CSidx.ih])))
% (6) greenhouse LTG range
fprintf("The LTG ranges between %.0f and %.0fC during greenhouse intervals\n",...
    min(ltgsat([CSidx.gh;CSidx.hh])),max(ltgsat([CSidx.gh;CSidx.hh])))
% (7) LTG range over the last 90 Ma
fprintf("Over the last 90 Ma, the LTG ranges between %.0f and %.0fC\n",...
    min(ltgsat(GTS.Average<90)),max(ltgsat(GTS.Average<90)))
% (8) SST LTG range over the last 90 Ma
fprintf("Over the last 90 Ma, the SST LTG ranges between %.0f and %.0fC\n",...
    min(ltgsst(GTS.Average<90)),max(ltgsst(GTS.Average<90)))
% (9) Average tropical temperatures during icehouse intervals
fprintf("The average tropical SAT ranges between %.0f and %.0fC (median = %.0f) during icehouse intervals\n",...
    min(tropsat([CSidx.ch;CSidx.ih])),max(tropsat([CSidx.ch;CSidx.ih])),...
    median(tropsat([CSidx.ch;CSidx.ih])))
% (10) Average tropical temperatures during greenhouse intervals
fprintf("The average tropical SAT ranges between %.0f and %.0fC (median = %.0f) during greenhouse intervals\n",...
    min(tropsat([CSidx.gh;CSidx.hh])),max(tropsat([CSidx.gh;CSidx.hh])),...
    median(tropsat([CSidx.gh;CSidx.hh])))
% (11) Tropical SSTs during the Turonian and PETM
fprintf("The average tropical SST during the Turonian was %.0f\n",...
    median(tropsst(GTS.Stage=="Turonian")))
fprintf("The average tropical SST during the PETM was %.0f\n",...
    median(tropsst(GTS.Stage=="PETM")))

% Phanerozoic GMST and atmospheric CO2
% (1) Phanerozoic correlation between CO2 and GMST
fprintf("The Phanerozoic correlation coefficient between GMST and CO2 is: %.2f\n",rA)
% (2) Cenozoic correlation between CO2 and GMST
fprintf("The Cenozoic correlation coefficient between GMST and CO2 is: %.2f\n",rC)
% (3) Paleozoic correlation between CO2 and GMST
fprintf("The Paleozoic correlation coefficient between GMST and CO2 is: %.2f\n",rP)
% (4) Mesozoic correlation between CO2 and GMST
fprintf("The Mesozoic correlation coefficient between GMST and CO2 is: %.2f (p-value = %.04f)\n",...
    rM,pvalM)
% (5) Median CO2 during coldhouse intervals
fprintf("The median CO2 during coldhouse intervals is: %.0f\n",...
    median(PhanerozoicCO2(intersect(CSidx.ih,allnocopse),:),[2,1]))
% (6) Median CO2 during transitional intervals
fprintf("The median CO2 during transitional intervals is: %.0f\n",...
    median(PhanerozoicCO2(intersect(CSidx.tr,allnocopse),:),[2,1]))
% (6) Median CO2 during hothouse intervals
fprintf("The median CO2 during hothouse intervals is: %.0f\n",...
    median(PhanerozoicCO2(intersect(CSidx.hh,allnocopse),:),[2,1]))
% (6) Median CO2 during hothouse intervals, excluding Mesozoic
fprintf("The median CO2 during hothouse intervals (excluding the Mesozoic) is: %.0f\n",...
    median(PhanerozoicCO2(setdiff(intersect(CSidx.hh,allnocopse),mesozoic),:),[2,1]))
% (7) Correlation coefficient between DFsolar+co2 and GMST
fprintf("The correlation coefficient between dFco2+sol and GMST is: %.2f\n",...
    corr(median(dFco2(allnocopse,:),2)+dFsun(allnocopse,:),...
    pgmst(allnocopse,p==50)))
% (8) Correlation coefficient between DFsolar+co2+surface and GMST
fprintf("The correlation coefficient between dFtot and GMST is: %.2f\n",...
    corr(median(dFco2(allnocopse,:),2)+dFsun(allnocopse,:)+...
    dFsa(allnocopse,:),pgmst(allnocopse,p==50)))
% (9) Phanerozoic ACS
fprintf('The Phanerozoic ACS is %.1f %s %.1f\n', mA,char(177),smA);
% (10) Cenozoic ACS
fprintf('The Cenozoic ACS is %.1f %s %.1f\n', mC,char(177),smC);
% (11) Paleozoic ACS
fprintf('The Paleozoic ACS is %.1f %s %.1f\n', mP,char(177),smP);
% (12) Paleozoic ACS
fprintf('The Cenozoic & Paleozoic ACS is %.1f %s %.1f\n', mCP,char(177),smCP);

