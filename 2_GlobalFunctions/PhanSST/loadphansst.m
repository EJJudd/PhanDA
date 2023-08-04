function Data = loadphansst(datafilename)

% Specify which fields are strings vs. numeric values in data file
stringfields = {'SampleID','SiteName','SiteHole','Formation','Country',...
        'ContinentOcean','Period','Stage','StagePosition','Biozone',...
        'ProxyType','ValueType','Taxon1','Taxon2','Taxon3','Environment',...
        'Ecology','CL','LeadAuthor','PublicationDOI','DataDOI'};
doublefields = {'MBSF','MCD','SampleDepth','ModLat','ModLon','Age',...
        'AgeFlag','ProxyValue','DiagenesisFlag','Mn','Fe','Sr','Mg','Ca',...
        'Cawtp','MgCa','SrCa','MnSr','NBS120c','Durango','MaximumCAI',...
        'ModWaterDepth','CleaningMethod','GDGT0','GDGT1','GDGT2','GDGT3',...
        'Cren','Crenisomer','BIT','dRI','MI','Year'};

% Save options
opts = detectImportOptions(datafilename);
opts = setvartype(opts,stringfields,'char');
opts = setvartype(opts,doublefields,'double');
opts = setvaropts(opts,stringfields,'FillValue','');

% Load data
Data = readtable(datafilename,opts);

for ii = 1:numel(stringfields)
    Data.(stringfields{ii}) = string(Data.(stringfields{ii}));
end

end