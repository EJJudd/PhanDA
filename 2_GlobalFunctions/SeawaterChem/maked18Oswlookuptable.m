function SeawaterLookup = maked18Oswlookuptable(Data, Preferences, ...
    AssumptionFiles, Nens, ExpNo, Ens, EnsMat, EnsMeta, Lat, Lon, ...
    GTS, filename)

% (1) Preload necessary information
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
SeawaterLookup = table(Data.Stage(idx),Data.PaleoLat(idx), ...
    Data.PaleoLon(idx), 'VariableNames', varNames);
SeawaterLookup = unique(SeawaterLookup,'rows');
SeawaterLookup = table2struct(SeawaterLookup,"ToScalar",true);
SeawaterLookup.d18Oswlocal = NaN(numel(SeawaterLookup.Stage),Nens);
SeawaterLookup.EnsembleIndex = NaN(numel(SeawaterLookup.Stage),Nens);
Nsites = numel(SeawaterLookup.Stage);

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
SeawaterLookup.tos = NaN(Nsites, Nens); 
SeawaterLookup.so = NaN(Nsites, Nens);

% (4) Cycle through the stages to estimate local d18Osw
stages = SeawaterLookup.Stage;
stages(contains(stages,'/')) = extractBefore(stages(contains(stages,'/')),'/');
s = unique(stages);
for ii = 1:numel(s)

    % (a) index the stage, the sites of that stage, and ensemble members
    sidx = find(ExpNo.Stage == s(ii));
    didx = find(stages == s(ii));
    expidx = [ExpNo.ExpNo1(sidx), ...
              ExpNo.ExpNo2(sidx)];
    ensidx = any(expno == expidx,2);
    SeawaterLookup.EnsembleIndex(didx,:) = ...
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
             [SeawaterLookup.PaleoLat(didx(jj)), ...
             SeawaterLookup.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rsoanom(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "soanom", ...
             [SeawaterLookup.PaleoLat(didx(jj)), ...
             SeawaterLookup.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rtosanom(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "tosanom", ...
             [SeawaterLookup.PaleoLat(didx(jj)), ...
             SeawaterLookup.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rdist(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "dist", ...
             [SeawaterLookup.PaleoLat(didx(jj)), ...
             SeawaterLookup.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rtos(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "tos", ...
             [SeawaterLookup.PaleoLat(didx(jj)), ...
             SeawaterLookup.PaleoLon(didx(jj))], ...
             'exclude', exclude));
        Rso(didx(jj),1) = min(ensusemeta.closestLatLon(...
             "so", ...
             [SeawaterLookup.PaleoLat(didx(jj)), ...
             SeawaterLookup.PaleoLon(didx(jj))], ...
             'exclude', exclude));

    end
    % (c) Extract model prior values from sites
    x.precipln_anomaly(didx,:) = ensusemat(Rpranom(didx),:);
    x.sss_anomaly(didx,:) = ensusemat(Rsoanom(didx),:);
    x.sst_anomaly(didx,:) = ensusemat(Rtosanom(didx),:);
    x.dist(didx,:) = ensusemat(Rdist(didx),:);
    x.lat(didx,:) = repmat(Lat(mod(Rtosanom(didx),ensusemeta.lengths(1))),1, Nens);
    SeawaterLookup.tos(didx,:) = ensusemat(Rtos(didx),:);
    SeawaterLookup.so(didx,:) = ensusemat(Rso(didx),:);
    SeawaterLookup.lat(didx,:) = repmat(Lat(mod(Rtosanom(didx),ensusemeta.lengths(1))),1, Nens);
    SeawaterLookup.lon(didx,:) = repmat(Lon(mod(Rtosanom(didx),ensusemeta.lengths(1))),1, Nens);
    fprintf('Completed Row Assignments for Stage %d/%d\n', ii, numel(s))

end

% % Fill in NaN so & soanom values in Scotese 06 runs
% [r, ~] = find(isnan(x.sss_anomaly));
% fillval1 = nanmedian(SeawaterLookup.so(r,:),2);
% fillval2 = nanmedian(x.sss_anomaly(r,:),2);
% idx = find(isnan(x.sss_anomaly));
% SeawaterLookup.so(idx) = fillval1;
% x.sss_anomaly(idx) = fillval2;

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
SeawaterLookup.d18Oswlocal = reshape(d18Oswlocal, Nsites, Nens);
    

% (7) Save file
if exist('filename')
    save(filename,'SeawaterLookup');
end

end
    
    