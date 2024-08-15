%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAT & SST LTG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 1: LOAD DATA
% Directory details
assdate = '21May2024';
assdir = ['/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/AssimilationOutputs/PhanerozoicDA_',assdate];
% Load data
% PART 1: LOAD DATA
load([assdir,'/OutputWorkspaces/','Output.mat'],"GMST","LTG","ItName","LTGsst")
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
LTGsst = combineruns(LTGsst,idx,2);

% Revise LTG, GMST, & GTS to account for combined stages
endsize = size(GMST,1) - numel(Preferences.combstages(1,:));
GTS = combinestages(GTS,"GTS",Preferences,endsize);
GMST = combinestages(GMST,"GMST",Preferences,endsize,1);
LTG = combinestages(LTG,"LTG",Preferences,endsize,2);
LTGsst = combinestages(LTGsst,"LTGsst",Preferences,endsize,2);

% Define the regions
poles = abs(Lat)>=66.5;
tropics = abs(Lat)<=23.5;

% Calculate percentiles
p = [16,50,84];
pgmst = cell2mat(cellfun(@(x) prctile(x,p), GMST, 'UniformOutput', false));

% subdivide median GMST by climate state quantile
q = prctile(pgmst(:,p==50),[20:20:100]);
CSidx.ih = find(pgmst(:,p==50)<=q(1));
CSidx.ch = find(pgmst(:,p==50)>q(1) & pgmst(:,p==50)<=q(2));
CSidx.tr = find(pgmst(:,p==50)>q(2) & pgmst(:,p==50)<=q(3));
CSidx.gh = find(pgmst(:,p==50)>q(3) & pgmst(:,p==50)<=q(4));
CSidx.hh = find(pgmst(:,p==50)>q(4));

%% Calculate the polar/tropical temps + gradients
ltgsat = cell(size(GMST));
ltgsst = cell(size(GMST));
tropsat = cell(size(GMST));
tropsst = cell(size(GMST));
polesat = cell(size(GMST));
polesst = cell(size(GMST));
meandim = 1;
for ii = 1:numel(GMST)
    xt = LTG{ii}; xt(~tropics,:) = NaN;
    xp = LTG{ii}; xp(~poles,:) = NaN;
    ltgsat{ii} = (latweightgmst(xt,meandim)-latweightgmst(xp,meandim))';
    tropsat{ii} = (latweightgmst(xt,meandim))';
    polesat{ii}= (latweightgmst(xp,meandim))';
    xt = LTGsst{ii}; xt(~tropics,:) = NaN;
    xp = LTGsst{ii}; xp(~poles,:) = NaN;
    ltgsst{ii} = (latweightgmst(xt,meandim)-latweightgmst(xp,meandim))';
    tropsst{ii} = (latweightgmst(xt,meandim))';
    polesst{ii} = (latweightgmst(xp,meandim))';
end

% Determine the temps/gradients by climate state
% Approach #1 - take the median LTG from each stage
LTGSAT.md.ih = cellfun(@(idx) median(ltgsat{idx}), num2cell(CSidx.ih));
LTGSAT.md.ch = cellfun(@(idx) median(ltgsat{idx}), num2cell(CSidx.ch));
LTGSAT.md.tr = cellfun(@(idx) median(ltgsat{idx}), num2cell(CSidx.tr));
LTGSAT.md.wh = cellfun(@(idx) median(ltgsat{idx}), num2cell(CSidx.gh));
LTGSAT.md.hh = cellfun(@(idx) median(ltgsat{idx}), num2cell(CSidx.hh));
TROPSAT.md.ih = cellfun(@(idx) median(tropsat{idx}), num2cell(CSidx.ih));
TROPSAT.md.ch = cellfun(@(idx) median(tropsat{idx}), num2cell(CSidx.ch));
TROPSAT.md.tr = cellfun(@(idx) median(tropsat{idx}), num2cell(CSidx.tr));
TROPSAT.md.wh = cellfun(@(idx) median(tropsat{idx}), num2cell(CSidx.gh));
TROPSAT.md.hh = cellfun(@(idx) median(tropsat{idx}), num2cell(CSidx.hh));
POLESAT.md.ih = cellfun(@(idx) median(polesat{idx}), num2cell(CSidx.ih));
POLESAT.md.ch = cellfun(@(idx) median(polesat{idx}), num2cell(CSidx.ch));
POLESAT.md.tr = cellfun(@(idx) median(polesat{idx}), num2cell(CSidx.tr));
POLESAT.md.wh = cellfun(@(idx) median(polesat{idx}), num2cell(CSidx.gh));
POLESAT.md.hh = cellfun(@(idx) median(polesat{idx}), num2cell(CSidx.hh));
LTGSST.md.ih = cellfun(@(idx) median(ltgsst{idx}), num2cell(CSidx.ih));
LTGSST.md.ch = cellfun(@(idx) median(ltgsst{idx}), num2cell(CSidx.ch));
LTGSST.md.tr = cellfun(@(idx) median(ltgsst{idx}), num2cell(CSidx.tr));
LTGSST.md.wh = cellfun(@(idx) median(ltgsst{idx}), num2cell(CSidx.gh));
LTGSST.md.hh = cellfun(@(idx) median(ltgsst{idx}), num2cell(CSidx.hh));
TROPSST.md.ih = cellfun(@(idx) median(tropsst{idx}), num2cell(CSidx.ih));
TROPSST.md.ch = cellfun(@(idx) median(tropsst{idx}), num2cell(CSidx.ch));
TROPSST.md.tr = cellfun(@(idx) median(tropsst{idx}), num2cell(CSidx.tr));
TROPSST.md.wh = cellfun(@(idx) median(tropsst{idx}), num2cell(CSidx.gh));
TROPSST.md.hh = cellfun(@(idx) median(tropsst{idx}), num2cell(CSidx.hh));
POLESST.md.ih = cellfun(@(idx) median(polesst{idx}), num2cell(CSidx.ih));
POLESST.md.ch = cellfun(@(idx) median(polesst{idx}), num2cell(CSidx.ch));
POLESST.md.tr = cellfun(@(idx) median(polesst{idx}), num2cell(CSidx.tr));
POLESST.md.wh = cellfun(@(idx) median(polesst{idx}), num2cell(CSidx.gh));
POLESST.md.hh = cellfun(@(idx) median(polesst{idx}), num2cell(CSidx.hh));


% Approach #2 - use the whole ensemble
LTGSAT.en.ih = cat(1,ltgsat{CSidx.ih});
LTGSAT.en.ch = cat(1,ltgsat{CSidx.ch});
LTGSAT.en.tr = cat(1,ltgsat{CSidx.tr});
LTGSAT.en.wh = cat(1,ltgsat{CSidx.gh});
LTGSAT.en.hh = cat(1,ltgsat{CSidx.hh});
TROPSAT.en.ih = cat(1,ltgsat{CSidx.ih});
TROPSAT.en.ch = cat(1,ltgsat{CSidx.ch});
TROPSAT.en.tr = cat(1,ltgsat{CSidx.tr});
TROPSAT.en.wh = cat(1,ltgsat{CSidx.gh});
TROPSAT.en.hh = cat(1,tropsat{CSidx.hh});
POLESAT.en.ih = cat(1,ltgsat{CSidx.ih});
POLESAT.en.ch = cat(1,ltgsat{CSidx.ch});
POLESAT.en.tr = cat(1,ltgsat{CSidx.tr});
POLESAT.en.wh = cat(1,ltgsat{CSidx.gh});
POLESAT.en.hh = cat(1,polesat{CSidx.hh});
LTGSST.en.ih = cat(1,ltgsst{CSidx.ih});
LTGSST.en.ch = cat(1,ltgsst{CSidx.ch});
LTGSST.en.tr = cat(1,ltgsst{CSidx.tr});
LTGSST.en.wh = cat(1,ltgsst{CSidx.gh});
LTGSST.en.hh = cat(1,ltgsst{CSidx.hh});
TROPSST.en.ih = cat(1,ltgsst{CSidx.ih});
TROPSST.en.ch = cat(1,ltgsst{CSidx.ch});
TROPSST.en.tr = cat(1,ltgsst{CSidx.tr});
TROPSST.en.wh = cat(1,ltgsst{CSidx.gh});
TROPSST.en.hh = cat(1,tropsst{CSidx.hh});
POLESST.en.ih = cat(1,ltgsst{CSidx.ih});
POLESST.en.ch = cat(1,ltgsst{CSidx.ch});
POLESST.en.tr = cat(1,ltgsst{CSidx.tr});
POLESST.en.wh = cat(1,ltgsst{CSidx.gh});
POLESST.en.hh = cat(1,polesst{CSidx.hh});

% Make gradient tables
mdtab = table('Size',[5,5],'VariableTypes',["string","double","string","double","string"],...
    'VariableNames',["State","SAT Median","SAT Range","SST Median","SST Range"]);
mdtab.State = ["Coldhouse","Coolhouse","Transitional","Warmhouse","Hothouse"]';
entab = table('Size',[5,5],'VariableTypes',["string","double","string","double","string"],...
    'VariableNames',["State","SAT Median","SAT Std","SST Median","SST Std"]);
entab.State = ["Coldhouse","Coolhouse","Transitional","Warmhouse","Hothouse"]';
stateorder = ["ih","ch","tr","wh","hh"];
for ii = 1:numel(stateorder)
    mdtab.("SAT Median")(ii) = round(median(LTGSAT.md.(stateorder{ii})));
    mdtab.("SAT Range")(ii) = [num2str(round(min(LTGSAT.md.(stateorder{ii})))),' - '...
        num2str(round(max(LTGSAT.md.(stateorder{ii}))))];
    mdtab.("SST Median")(ii) = round(median(LTGSST.md.(stateorder{ii})));
    mdtab.("SST Range")(ii) = [num2str(round(min(LTGSST.md.(stateorder{ii})))),' - '...
        num2str(round(max(LTGSST.md.(stateorder{ii}))))];
    entab.("SAT Median")(ii) = round(median(LTGSAT.en.(stateorder{ii})));
    entab.("SAT Std")(ii) = [num2str(round(prctile(LTGSAT.en.(stateorder{ii}),16))),' - '...
        num2str(round(prctile(LTGSAT.en.(stateorder{ii}),84)))];
    entab.("SST Median")(ii) = round(median(LTGSST.en.(stateorder{ii})));
    entab.("SST Std")(ii) = [num2str(round(prctile(LTGSST.en.(stateorder{ii}),16))),' - '...
        num2str(round(prctile(LTGSST.en.(stateorder{ii}),84)))];
end





