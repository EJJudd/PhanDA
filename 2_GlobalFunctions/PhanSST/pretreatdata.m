function Data = pretreatdata(Data, Preferences, Corrections, ExpNo)

% Most preferences are specified in the parse data stage, while most
% corrections are applied to the Ye values, but for logistal reasons,
% the standardization of phosphate values takes place on the data side, as
% a pretreatment. This step also includes other pretreatment options, like
% paritioning between foram & macrofossil carbonate data opting to parse
% data by nearest grid cell vs rounded coords vs raw coords

load("GTS2020_PETM.mat")

% Pretreatment 0: Remove Cambrian data
Data(Data.Period == "Cambrian",:) = [];

% Pretreatment 1: NBS120c Correction
idx = find(strcmpi(Data.ProxyType, 'd18p'));
Data.ProxyValue(idx) = Data.ProxyValue(idx) - ...
        (Data.NBS120c(idx) - Corrections.NBS120c);
    
% Pretreatment 2: Durango Correction
idx = find(strcmpi(Data.ProxyType, 'd18p') & ...
    strcmpi(Data.AnalyticalTechnique, 'SIMS'));
Data.DiagenesisFlag(idx(isnan(Data.Durango(idx)))) = 1;
Data.ProxyValue(idx) = Data.ProxyValue(idx) - ...
        (Data.Durango(idx) - Corrections.Durango);

% Pretreatment 3: SIMS correction
Data.ProxyValue(idx) = Data.ProxyValue(idx) - Corrections.SIMS;

% Pretreatment 4: Parse foram vs macrofossil data
Data.ProxyType(strcmpi(Data.ProxyType,'d18c') & strcmpi(Data.Taxon1, 'pf')) = {'d18cforam'};
Data.ProxyType(strcmpi(Data.ProxyType,'d18c')) = {'d18cmacro'};

% Pretreatment 5: Indicate PETM data
Data = addpetm(Data, Preferences);

% Pretreatment 6: Assign paleo-coordinates
load("PaleoCoordinates.mat")
Data = expno2paleolat(Data, PaleoCoordinates, ExpNo);

% Pretreatment 7: Adjust lat/lon based on method
if strcmpi(Preferences.coordmethod, "gridcell")
    load("HadCM3Coordinates.mat", "Lat", "Lon")
    Data.PaleoLat = Lat(dsearchn(Lat, Data.PaleoLat));
    Data.PaleoLon = Lon(dsearchn(Lon, Data.PaleoLon));  
elseif strcmpi(Preferences.coordmethod, "round")
    Data.PaleoLat = round(Data.PaleoLat);
    Data.PaleoLon = round(Data.PaleoLon);
elseif strcmpi(Preferences.coordmethod, "raw")
else
    warning("Invalid coord method preference provided")
end

% Pretreatment 8: Add paleodepths
Data = addpaleodepths(Data,ExpNo);

% Pretreatment 9: Update stages
for ii = 1:size(Preferences.combstages,2)
    Data.Stage(Data.Stage == GTS.Stage(Preferences.combstages(1,ii))) = ...
        strcat(GTS.Stage(Preferences.combstages(1,ii)),"/",...
        GTS.Stage(Preferences.combstages(2,ii)));
    Data.Stage(Data.Stage == GTS.Stage(Preferences.combstages(2,ii))) = ...
        strcat(GTS.Stage(Preferences.combstages(2,ii)),"/",...
        GTS.Stage(Preferences.combstages(1,ii)));
end

end