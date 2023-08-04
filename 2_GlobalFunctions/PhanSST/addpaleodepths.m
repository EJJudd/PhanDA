function Data = addpaleodepths(Data,ExpNo)

PaleoWaterDepth = NaN(height(Data),1);
Data = addvars(Data,PaleoWaterDepth,'After','ModWaterDepth');

bathdir = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs/ModelOutputs/Bathymetry/';
bathfiles = strcat(bathdir,string(extractfield( dir( [bathdir,'*.nc']), 'name'))');

load("HadCM3Coordinates.mat");
[Lon,Lat] = meshgrid(Lon,Lat);

idx = find(Data.ProxyType == "mg");

uniqueagecoor = unique( table(Data.Stage(idx), Data.PaleoLat(idx), Data.PaleoLon(idx), ...
    'VariableNames', {'Stage','PaleoLat','PaleoLon'}), 'rows');

s = string(unique(uniqueagecoor.Stage));
for ii = 1:numel(s)
    sidx = find(string(ExpNo.Stage) == s(ii));
    uidx = find(string(uniqueagecoor.Stage) == s(ii));
    expnos = [ExpNo.ExpNo1(sidx);ExpNo.ExpNo2(sidx)];
    bath1 = ncread(bathfiles(expnos(1)),'depth');
    bath2 = ncread(bathfiles(expnos(2)),'depth');
    for jj = 1:numel(uidx)
        depth = nanmean([bath1(Lat == uniqueagecoor.PaleoLat(uidx(jj)) & ...
                              Lon == uniqueagecoor.PaleoLon(uidx(jj))), ...
                         bath2(Lat == uniqueagecoor.PaleoLat(uidx(jj)) & ...
                              Lon == uniqueagecoor.PaleoLon(uidx(jj)))]);
        depth(isnan(depth)) = 0;
        Data.PaleoWaterDepth(Data.ProxyType == "mg" & Data.Stage == s(ii) & ...
            Data.PaleoLat == uniqueagecoor.PaleoLat(uidx(jj)) & ...
            Data.PaleoLon == uniqueagecoor.PaleoLon(uidx(jj))) = depth;
    end
end

end
