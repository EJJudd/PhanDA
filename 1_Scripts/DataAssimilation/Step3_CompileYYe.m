%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSIMILATION STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Compile Seawater Lookup table, Y, & Ye  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 05/21 (E. Judd)
% Last updated: 08/22 (E. Judd)

% Notes: 

%   -->  This script is compatable with the most recent Dash release
%        (https://github.com/JonKing93/DASH/releases/tag/v4.2.2)

%   -->  This script is divided into x Parts
%        + Part 0: Define (and make) the assimilation directories
%        + Part 1: Load data
%        + Part 2: Pre-treat data
%        + Part 3: Make (or append, load existing) seawater lookup table
%        + Part 4: Parse data by stage and calculate Y, Ye
%   -->  Variable naming Convention:
%        + Permanant variables: begin with uppercase, camel case throughout
%        + Temporary or overwriten variables: all lowercase


%% PART 0: DEFINE (AND MAKE) THE ASSIMILATION DIRECTORIES

% ----- UPDATE TO MATCH YOU FILEPATHS AND FOLDER CONFIGURATION -----
% Root directories (in my case, I store the input data -- i.e., the proxy 
% data and model priors -- on iCloud, and save the assimilation outputs on
% OneDrive)
iCloud = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs';
OneDrive = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC';

% Model Prior Ensemble:
% Ensemble Directory (stored on iCloud)
EnsDir = [iCloud,'/ModelOutputs/AssimilationFiles_beta'];
% Ensemble filename
enfilename = [EnsDir,'/HadCM3_all.ens'];

% PhanSST Data:
% Data Directory (stored on OneDrive)
DataDir = [iCloud,'/PhanSST_RawData/OutputFiles'];
% data filename
datafilename = [DataDir,'/PhanSST_v001_06Jun2023.csv'];

% Assimilaion Directory:
% (where input, output, and figure files are saved)
% I name these directories by date for easy archiving
ad = strrep(date,'-','');
AssDir = [OneDrive,'/AssimilationOutputs/PhanerozoicDA_',ad];
mkdir(AssDir)
% Define and create the assimilation directory subfolders
InputDir = [AssDir,'/InputWorkspaces'];
OutputDir = [AssDir,'/OutputWorkspaces'];
FigDir = [AssDir,'/OutputFigs'];
mkdir(InputDir)
mkdir(OutputDir)
mkdir(FigDir)
% ----- DO NOT UPDATE PAST THIS LINE -----


%% PART 1: LOAD DATA

% (b) Read in Files
% PhanSST data file
Data = loadphansst(datafilename);

% Prior Ensemble
Ens = ensemble(enfilename);
EnsMat = Ens.load;
EnsMeta = Ens.metadata;
% Geologic time scale, 2020
load("GTS2020_PETM.mat","GTS");
% Experiment number/plate rotation age file
load("ExpNo.mat", "ExpNo")
% Latitude & longitude vectors
load("HadCM3Coordinates.mat","Lat","Lon")


%% PART 2: PRE-TREAT DATA
% (a) Define DA range by stage number (Holcene = 1, Tremadocian = 91)
DArange = 1:91;

% (b) Decision making regarding pretrement of data and calculation of Ye values          
Preferences.coordmethod = "gridcell";
Preferences.agemethod = "stage";
Preferences.diagenesisflag = "exclude";
Preferences.maxCAI = 4;
Preferences.layer = "surfaceonly";
Preferences.environment = ["estuarine","freshwater","tidal flat","lagoon","restricted"];
Preferences.combstages = [22,45,51,56,75,83;...
                          23,46,52,57,76,84];
                      Preferences.PETM = "body";
Preferences.d18a.eq = "Grossman1986"; Preferences.d18a.min = "aragonite";
Preferences.d18p.eq = "Lecuyer2013"; Preferences.d18p.min = "phosphate";
Preferences.d18c.eq.nobel = "Kim1997"; Preferences.d18c.min = "calcite";
Preferences.d18c.eq.bel = "Daeron2019";

% (c) Decision making regarding corrections for Ye values
    % NOTE: if using the Daeron 2019 equation for the belemnite data, set
    % the belemnite correction to 0; otherwise, it is recommeneded to use a
    % correction of 1.5;
Corrections.globalsw = "theory";
Corrections.snowballearth = true;
Corrections.phcorrection = true;
Corrections.NBS120c = 21.7;
Corrections.Durango = 9.8;
Corrections.SIMS = .5;
Corrections.belemnite = 0;

% (d) pre-treat data
% Manually adjust modern coordinates of the Iranian data from Angiolini et 
% al. (2009) away from plate boundary prior to rotating
Data.ModLat(Data.PublicationDOI=="https://doi.org/10.1144/0016-76492008-096R" & ...
    Data.Country == "Iran") = 36.01;
Data.ModLon(Data.PublicationDOI=="https://doi.org/10.1144/0016-76492008-096R" & ...
    Data.Country == "Iran") = 51.48;

% NOTE: Most preferences are specified in the parse data stage, while most
% corrections are applied to the Ye values, but for logistal reasons,
% the standardization of phosphate values takes place on the data side, as
% a pretreatment. This step also includes other pretreatment options, 
% including:
% + seperating out the PETM data
% + adding bathymetry values to Mg/Ca data
% + adding paleocoordinats
% + paritioning between foram & macrofossil carbonate data
% + option to parse data by nearest grid cell or rounded coords
% + updating stage names to reflect combined stages
Data = pretreatdata(Data, Preferences, Corrections, ExpNo);

% (e) Manual data adjustments
Data = removedflaggedsites(Data);

% (f) Load assumption files
AssumptionFiles = loadassumptionfiles(Corrections);

% (g) Save
save([InputDir,'/Data.mat'],'Data')

%%  PART 3: MAKE LOCAL SEAWATER LOOKUP TABLE

% Name the file
% Note: Saving the output to the GlobalFiles folder makes it globally
%       accessible
filename = [OneDrive,'/Code/DataAssimilation/4_GlobalFiles/PSMs/seawaterchem/SeawaterLookup.mat'];

% OPTION A: MAKE A NEW TABLE
% Syntax: (uncomment to use)
%     % How many ensemble members per scenario?
%     Nens = 80;
%     % Make lookup table
%     SeawaterLookup = maked18Oswlookuptable(Data, Preferences, ...
%         AssumptionFiles, Nens, ExpNo, Ens, EnsMat, EnsMeta, Lat, Lon, ...
%         GTS, filename);


% OPTION B: append existing table
% Note: if archive = true, then the old and new filenames can be the same;
%       the old lookup table will be moved to the archive folder and
%       renamed to include the data of the move.
% Syntax: (uncomment to use)
%     Nens = 80;
%     filenameold = filename;
%     filenamenew = filename;
%     archive = true;
%     archivedir = [OneDrive,'/Code/DataAssimilation/Archive/4_GlobalFiles'];
%     SeawaterLookup = appendd18Oswlookuptable(Data, Preferences, ...
%         AssumptionFiles, Nens, ExpNo, Ens, EnsMat, EnsMeta, ...
%         filenameold, filenamenew, archive, archivedir);

% OPTION C: load existing table
% Syntax: (uncomment to use)
%     % If not saved as a global file:
%     load(filename)
%     % If saved as a global file:
    load("SeawaterLookup.mat", "SeawaterLookup")

%%  PART 4: MAKE GLOBAL pH LOOKUP TABLE
filename = [OneDrive,'/Code/DataAssimilation/4_GlobalFiles/PSMs/seawaterchem/pHLookup.mat'];
% OPTION A: MAKE A NEW TABLE
% Syntax: (uncomment to use)
%    % How many ensemble members per scenario?
%    Nens = 80;    
%    % Make lookup table
%    pHLookup = makepHlookuptable(DArange, Nens, ExpNo, EnsMat, EnsMeta, ...
%         GTS, Lat, Lon, filename);
% OPTION B: load existing table
% Syntax: (uncomment to use)
%     % If not saved as a global file:
%     load(filename)
%     % If saved as a global file:
     load("pHLookup.mat", "pHLookup")

%% PART 5: PARSE DATA & CALCULATE Y/YE

% How should the ensemble members be indexed?
dims = string(EnsMeta.ensembleDimensions{1});
expno = EnsMeta.ensemble{1,1}{1, dims == "expno"};
% What labels should be used for the data, Y, and Ye structures?
stagelab = strrep( strcat( ...
        'S', num2str(DArange','%02.f'), '_', GTS.Stage(DArange)), ' ', '');
% Do you want to print the progress of each Y, Ye calculation?
yeverbose = true;
% Do you want to print the progress after each stage?
stageverbose = true;


% (3)Assemble unique proxy data (UPD), Y, Ye, and assumption values
% To assemble the Y/Ye Files:

for a = DArange
%       (a) Parse ensemble to use
    sidx = find(ExpNo.Stage == GTS.Stage(a));
    expidx = [ExpNo.ExpNo1(sidx),ExpNo.ExpNo2(sidx)];
    ensidx = any(expno == expidx,2);
    ensuse = Ens.useMembers(ensidx);
    ensusemat = EnsMat(:,ensidx);
    ensusemeta = ensuse.metadata;

%       (b) Parse data by time interval and proxy
    [~,UPD.(stagelab{a})] = parsedata( ...
        Data, Preferences, GTS.Stage(a));

%       (d) Get inovation values and save assumption information
    [Y.(stagelab{a}), Ye.(stagelab{a}), Assumptions.(stagelab{a})] = ...
        estimateYYe(UPD.(stagelab{a}), ensuse, ensusemeta, a, ...
          GTS.Stage(a), SeawaterLookup, pHLookup, AssumptionFiles, ...
          Corrections, Preferences, GTS, yeverbose);

%       (e) Print progress (if option selected)        
    message = sprintf('Completed %s (Stage %d/%d)',GTS.Stage(a),a,numel(DArange));
    printprogress(message, stageverbose)

end

filename = [InputDir,'/YYe.mat'];
save(filename,'UPD','Y', 'Ye', 'Assumptions','Preferences','Corrections') 



