%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSIMILATION STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAKE GRIDFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 05/21 (E. Judd)
% Last updated: 08/23 (E. Judd)

% Notes: 
%   -->  This script is compatable with the most recent Dash release
%        (https://github.com/JonKing93/DASH/releases/tag/v4.2.2)
%   -->  This script builds a single gridfile that houses all HadCM3 models

clearvars;
close all;
clc;

%% MAKE GRIDFILE

% (1) Define (and/or make) directories & model files
% ----- UPDATE TO MATCH YOU FILEPATHS AND FOLDER CONFIGURATION -----
% Root directory (in my case, iCloud)
iCloud = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs';
% Model directory
mdldir = [iCloud,'/ModelOutputs/ExpFiles'];
% Gridfile directory (make if doesn't exist)
gridfiledir = [iCloud,'/ModelOutputs/AssimilationFiles'];
mkdir(gridfiledir)
% ----- DO NOT UPDATE PAST THIS LINE -----

% (2) Exctract model file names
allfiles = extractfield(dir([mdldir,'/*.nc']),'name')';
allfiles = strcat(mdldir, '/', allfiles);

% (3) Define grid metadata        
lon = ncread( string(allfiles(1)), 'lon' );
lat = ncread( string(allfiles(1)), 'lat' );
lev = [1];
time = [1:12]';
variable = ["tas";"tos";"so";"pr";"tosanom";"soanom";"prlnanom";"dist"];
expno = [1:109]';
exprun = ["scotesesolara";"scotese02";"scotese06";"scotese07";"scotese08";"scotesespinupa";"scotese2co2a";"scotese4co2a"];
expmean = [1:5]';
attributes = struct('Model', 'HadCM3', 'Grid', '2.5lat x 3.75lon');

% (4) Initialize grid file
metadata = gridMetadata( "lon", lon, "lat", lat, "lev", lev, ...
    "time", time, "var", variable, "expno", expno, "exprun", exprun, ...
     "expmean", expmean, "attributes", attributes);
overwrite = true;
HadCM3 = gridfile.new([gridfiledir,'/HadCM3_all.grid'], metadata, overwrite);


% (5) Add model files to the grid
% Asign the dimensional order for the files
% For standard variables:
dimorder = ["lat","lon","lev","time"];
% For anomaly variables
dimorderanom = ["lat","lon"];
% Loop through each simulation file to add data
for f = 1:numel( allfiles )

% Extract variable
expinfo = strsplit(string(allfiles(f)),"_");
exprun = sprintf("%s%s", "scotese", expinfo(2));
expno = str2double(expinfo(3));
expmean = str2num(strrep(expinfo(6),"M",""));
variable = regexprep(expinfo(7),".nc","");

% Add the data source, parsed by variable file
if variable == "pr"
    sourceMetadata = gridMetadata("lon",lon, "lat",lat, "lev", lev, ...
         "time", time, "exprun", exprun, "var", variable,...
          "expno", expno, "expmean", expmean);
    HadCM3.add( "nc", string(allfiles(f)), "pr", dimorder, sourceMetadata);
    sourceMetadata = gridMetadata("lon",lon, "lat",lat, "lev", lev, ...
         "time", 1, "exprun", exprun, "var", "prlnanom",...
          "expno", expno, "expmean", expmean);
    HadCM3.add( "nc", string(allfiles(f)), "prln_anomaly", dimorderanom, sourceMetadata);
elseif variable == "so"
    sourceMetadata = gridMetadata("lon",lon, "lat",lat, "lev", lev, ...
         "time", time, "exprun", exprun, "var", variable,...
          "expno", expno, "expmean", expmean);
    HadCM3.add( "nc", string(allfiles(f)), "so", dimorder, sourceMetadata);
    sourceMetadata = gridMetadata("lon",lon, "lat",lat, "lev", lev, ...
         "time", 1, "exprun", exprun, "var", "soanom",...
          "expno", expno, "expmean", expmean);
    HadCM3.add( "nc", string(allfiles(f)), "so_anomaly", dimorderanom, sourceMetadata);
elseif variable == "tas"
    sourceMetadata = gridMetadata("lon",lon, "lat",lat, "lev", lev, ...
         "time", time, "exprun", exprun, "var", variable,...
          "expno", expno, "expmean", expmean);
    HadCM3.add( "nc", string(allfiles(f)), "tas", dimorder, sourceMetadata);
elseif variable == "tos"
    sourceMetadata = gridMetadata("lon",lon, "lat",lat, "lev", lev, ...
         "time", time, "exprun", exprun, "var", variable,...
          "expno", expno, "expmean", expmean);
    HadCM3.add( "nc", string(allfiles(f)), "tos", dimorder, sourceMetadata);
    sourceMetadata = gridMetadata("lon",lon, "lat",lat, "lev", lev, ...
         "time", 1, "exprun", exprun, "var", "tosanom",...
          "expno", expno, "expmean", expmean);
    HadCM3.add( "nc", string(allfiles(f)), "tos_anomaly", dimorderanom, sourceMetadata);
    sourceMetadata = gridMetadata("lon",lon, "lat",lat, "lev", lev, ...
         "time", 1, "exprun", exprun, "var", "dist",...
          "expno", expno, "expmean", expmean);
    HadCM3.add( "nc", string(allfiles(f)), "distfromcoast", dimorderanom, sourceMetadata);
else
    fprintf("Error: no data added: %d\n",f)
end

if mod(f,250) == 0
    disp(f)
end

end
