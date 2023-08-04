function SeawaterLookup = appendd18Oswlookuptable(Data, Preferences, ...
    AssumptionFiles, Nens, ExpNo, Ens, EnsMat, EnsMeta, ...
    filenameold, filenamenew, archive, archivedir)

% (1) Preload necessary information
load(filenameold, "SeawaterLookup")
load("GTS2020_PETM.mat", "GTS")
load("HadCM3Coordinates.mat", "Lat", "Lon")
[Lat,Lon] = meshgrid(Lat,Lon);
Lat = Lat(:); Lon = Lon(:);
dims = string(EnsMeta.ensembleDimensions{1});
expno = EnsMeta.ensemble{1,1}{1, dims == "expno"};

% (2) Make structure with unique d18O lat, lon, and stages
% (a) Name vars
varNames = ["Stage","PaleoLat","PaleoLon"];
% (b) index data to include
if strcmpi(Preferences.diagenesisflag,'include')
    idx = find(contains(Data.ProxyType,'d18'));
elseif strcmpi(Preferences.diagenesisflag,'exclude')
    idx = find(contains(Data.ProxyType,'d18') & Data.DiagenesisFlag ~= 1);
end
% (c) make seperate tables for each ensemble scenario, find unique coors/
%     ages, convert tables to structures, & initialize seawater variable
dlu1 = table(SeawaterLookup.Stage, SeawaterLookup.PaleoLat, ...
    SeawaterLookup.PaleoLon, 'VariableNames', varNames);
dlu2 = table(Data.Stage(idx),Data.PaleoLat(idx), ...
    Data.PaleoLon(idx), 'VariableNames', varNames);
dlu2 = unique(dlu2,'rows');
dlu2 = setdiff(dlu2, dlu1);
dlu2 = table2struct(dlu2,"ToScalar",true);
dlu2.d18Oswlocal = NaN(numel(dlu2.Stage),Nens);
Nsites = numel(dlu2.Stage);

% (3) Find the appropriate variable values
%  Preallocate Matrices:
%  Row index vectors
Rpranom = NaN(Nsites, 1); 
Rsoanom = NaN(Nsites, 1);
Rtosanom = NaN(Nsites, 1); 
Rdist = NaN(Nsites, 1);
Rtos = NaN(Nsites, 1); 
Rso = NaN(Nsites, 1);
%  Prior values at proxy site matrices
x.lat = NaN(Nsites, Nens); 
x.precipln_anomaly = NaN(Nsites, Nens);
x.sss_anomaly = NaN(Nsites, Nens); 
x.sst_anomaly = NaN(Nsites, Nens);
x.dist = NaN(Nsites, Nens); 
dlu2.tos = NaN(Nsites, Nens); 
dlu2.so = NaN(Nsites, Nens);

% (4) Cycle through the stages to estimate local d18Osw
stages = dlu2.Stage;
stages(contains(stages,'/')) = extractBefore(stages(contains(stages,'/')),'/');
s = unique(stages);

for ii = 1:numel(s)
    % (a) index the stage, the sites of that stage, and ensemble members
    sidx = find(ExpNo.Stage == s(ii));
    didx = find(stages == s(ii));
    expidx = [ExpNo.ExpNo1(sidx), ...
              ExpNo.ExpNo2(sidx)];
    ensidx = any(expno == expidx,2);
    dlu2.EnsembleIndex(didx,:) = ...
        repmat(find(ensidx==1)',numel(didx),1);
    ensuse = Ens.useMembers(ensidx);
    ensusemat = EnsMat(:,ensidx);
    ensusemeta = EnsMeta.removeMembers(~ensidx);
    exclude = mean(load(ensuse.useVariables('tosanom')),2);
    exclude(~isnan(exclude)) = false;
    exclude(isnan(exclude)) = true;
    exclude = logical(exclude);

    % (b) Determine rows of variables used to estimate local d18sw
    for jj = 1:numel(didx)
        Rpranom(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "prlnanom", ...
             [dlu2.PaleoLat(didx(jj)), ...
             dlu2.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rsoanom(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "soanom", ...
             [dlu2.PaleoLat(didx(jj)), ...
             dlu2.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rtosanom(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "tosanom", ...
             [dlu2.PaleoLat(didx(jj)), ...
             dlu2.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rdist(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "dist", ...
             [dlu2.PaleoLat(didx(jj)), ...
             dlu2.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rtos(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "tos", ...
             [dlu2.PaleoLat(didx(jj)), ...
             dlu2.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rso(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "so", ...
             [dlu2.PaleoLat(didx(jj)), ...
             dlu2.PaleoLon(didx(jj))], ...
             'exclude', exclude));
    end
    % (c) Extract model prior values from sites
    x.precipln_anomaly(didx,:) = ensusemat(Rpranom(didx),:);
    x.sss_anomaly(didx,:) = ensusemat(Rsoanom(didx),:);
    x.sst_anomaly(didx,:) = ensusemat(Rtosanom(didx),:);
    x.dist(didx,:) = ensusemat(Rdist(didx),:);
    x.lat(didx,:) = repmat(Lat(mod(Rtosanom(didx),ensusemeta.lengths(1))),1, Nens);
    dlu2.tos(didx,:) = ensusemat(Rtos(didx),:);
    dlu2.so(didx,:) = ensusemat(Rso(didx),:);
    dlu2.lat(didx,:) = repmat(Lat(mod(Rtosanom(didx),ensusemeta.lengths(1))),1, Nens);
    dlu2.lon(didx,:) = repmat(Lon(mod(Rtosanom(didx),ensusemeta.lengths(1))),1, Nens);
    fprintf('Completed Row Assignments for Stage %d/%d (%s)\n', ii, numel(s),string(s(ii)))
end

% Fill in inf values (unclear why they appear?)
x.precipln_anomaly(isinf(x.precipln_anomaly)) = NaN;
[r, ~] = find(isnan(x.precipln_anomaly));
fillval = nanmedian(x.precipln_anomaly(r,:),2);
x.precipln_anomaly(isnan(x.precipln_anomaly)) = fillval; 

% (5) Estimate local seawater
x.lat = x.lat(:);
x.precipln_anomaly = x.precipln_anomaly(:);
x.sss_anomaly = x.sss_anomaly(:);
x.sst_anomaly = x.sst_anomaly(:);
x.dist = x.dist(:);
x = struct2table(x);
d18Oswlocal = predictd18Osw(x,AssumptionFiles,1);
dlu2.d18Oswlocal = reshape(d18Oswlocal, Nsites, Nens);

% (6) append new table to old table
load(filenameold)
    
SeawaterLookup.Stage = ...
    [SeawaterLookup.Stage; dlu2.Stage];
SeawaterLookup.PaleoLat = ...
    [SeawaterLookup.PaleoLat; dlu2.PaleoLat];
SeawaterLookup.PaleoLon = ...
    [SeawaterLookup.PaleoLon; dlu2.PaleoLon];
SeawaterLookup.d18Oswlocal = ...
    [SeawaterLookup.d18Oswlocal; dlu2.d18Oswlocal];
SeawaterLookup.EnsembleIndex = ...
    [SeawaterLookup.EnsembleIndex; dlu2.EnsembleIndex];
SeawaterLookup.tos = ...
    [SeawaterLookup.tos; dlu2.tos];
SeawaterLookup.so = ...
    [SeawaterLookup.so; dlu2.so];
SeawaterLookup.lat = ...
    [SeawaterLookup.lat; dlu2.lat];
SeawaterLookup.lon = ...
    [SeawaterLookup.lon; dlu2.lon];

% (7) Save new table
if archive == true
    archivename = strsplit(filenameold,'/');
    archivename = strcat(archivedir, '/', strrep(archivename{end},'.mat',''), ...
        '_', regexprep(date,'-',''), '.mat');
    movefile(filenameold,archivename)
end

if exist('filenamenew', 'var')
    save(filenamenew,'SeawaterLookup');
end


end
    
    