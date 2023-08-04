%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSIMILATION STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Compile Seawater Lookup table, Y, & Ye  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 05/21 (E. Judd)
% Last updated: 08/22 (E. Judd)

% Notes: 

%   -->  This script is compatable with the most recent Dash release
%        (https://github.com/JonKing93/DASH/releases/tag/v4.0.0-beta-2)

%   -->  This script is divided into x Parts
%        + Part 0: Define (and make) the assimilation directories
%        + Part 1: Load data
%        + Part 2: Pre-treat data
%        + Part 3: Make (or append, load existing) seawater lookup table
%        + Part 4: Parse data by stage and calculate Y, Ye
%   -->  Variable naming Convention:
%        + Permanant variables: begin with uppercase, camel case throughout
%        + Temporary or overwriten variables: all lowercase
%   -->  Graveyard code at end of this script to:
%        + Assign plate ages to the ExpNo file 
%   -->  adddist is a workaround to an error in haversine, it appears in
%        the following functions:
%        + haversine
%        Once this error is fixed, remove the adddist command
%



%% PART 0: DEFINE (AND MAKE) THE ASSIMILATION DIRECTORIES

% Home Directory (root for all other directories)
OneDrive = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC';
iCloud = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs';
% Ensemble Directory (stored on icloud)
EnsDir = [iCloud,'/ModelOutputs/AssimilationFiles_beta'];
DataDir = [iCloud,'/PhanSST_RawData/OutputFiles'];
% Assimilaion Directory (where input, output, and figure files are saved)
% Assimilation date:
% New assimilation:
 ad = strrep(date,'-','');
% Old assimilation
% ad = '17Jul2023';

AssDir = [OneDrive,'/AssimilationOutputs/PhanerozoicDA_',ad];
InputDir = [AssDir,'/InputWorkspaces'];
OutputDir = [AssDir,'/OutputWorkspaces'];
FigDir = [AssDir,'/OutputFigs'];

mkdir(AssDir)
mkdir(InputDir)
mkdir(OutputDir)
mkdir(FigDir)


%% PART 1: LOAD DATA

% (a) Filepaths
% Ensemble file
enfilename = [EnsDir,'/HadCM3_all.ens'];
datafilename = [DataDir,'/PhanSST_v001_06Jun2023.csv'];

% (b) Read in Files
% PhanSST data file
Data = loadphansst(datafilename);
% Remove Cambrian data
Data(Data.Period == "Cambrian",:) = [];

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
% (a) Define DA range (Holcene = 1)
DArange = 1:91;

% (b) Decision making regarding pretrement of data and calculation of Ye values          
Preferences.coordmethod = "gridcell";
Preferences.agemethod = "stage";
Preferences.diagenesisflag = "exclude";
Preferences.maxCAI = 4;
Preferences.layer = "surfaceonly";
Preferences.environment = ["estuarine","freshwater","tidal flat","lagoon"];
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
% Manually adjust data of Iranian data away from plate boundary
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
Data.DiagenesisFlag(Data.ProxyValue<-10) = 1; 
% --> (1) Gyanyima data
%         (from Garbelli et al., 2016)
Data.PaleoLat(Data.SiteName == "Gyanyima") = -40;
Data.PaleoLon(Data.SiteName == "Gyanyima") = 63.75;
% --> (2) Zanclean ODP978 data
%         (from Khelif et al., 2009 & 2014)
Data.DiagenesisFlag(Data.SiteName=="ODP978"&Data.Stage=="Zanclean") = 1;
% --> (3) Serpukhovian Askyn data 
%         (from Bruckschen et al., 2001 & Brand and Bruckschen, 2002)
Data.DiagenesisFlag(Data.SiteName=="Askyn"&Data.Stage=="Serpukhovian") = 1;
% --> (4) M'rirt data 
%         (from Le Houedec et al., 2013)
Data.DiagenesisFlag(Data.SiteName=="M'rirt") = 1;
% --> (5) d18p data from Xikou
%         (from Wang et al., 2020)
Data.DiagenesisFlag(Data.SiteName=="Xikou"&Data.ProxyType=="d18p",:) = 1;
% --> (5) data from the Kapp Starostin formation
%         (from Mii et al., 1997)
Data.DiagenesisFlag(Data.Formation=="Kapp Starostin") = 1;
% --> (6) d18cmacro deglacial/interglacial data
%         (from Frank et al., 2015)
Data.DiagenesisFlag(Data.PublicationDOI=="https://doi.org/10.1016/j.palaeo.2014.11.016" & ...
    (Data.Stage=="Sakmarian" | Data.Stage=="Artinskian") & Data.ProxyValue<-1.5) = 1;
% --> (7) d18cmacro data from Fayetteville Fm.
%         (from Grossman et al., 2008)
Data.PaleoLat(Data.Formation=="Fayetteville") = -12.5;
Data.PaleoLon(Data.Formation=="Fayetteville") = -33.75;
% --> (7) Serpukhovian d18p data from Tellego
%         (from Wang et al., 2020)
Data.PaleoLat(Data.SiteName=="Tellego"&Data.Stage=="Serpukhovian") = -12.5;
Data.PaleoLon(Data.SiteName=="Tellego"&Data.Stage=="Serpukhovian") = 7.5;
% --> (8) Lanadian, early Carnian, and early Norian
%         (from Hornung et al., 2007)
Data.DiagenesisFlag(contains(Data.Stage,"Ladinian")&Data.LeadAuthor=="Hornung") = 1;
Data.DiagenesisFlag(Data.Stage=="Carnian"&Data.StagePosition=="early"&Data.LeadAuthor=="Hornung") = 1;
Data.DiagenesisFlag(Data.Stage=="Norian"&Data.StagePosition=="early"&Data.LeadAuthor=="Hornung") = 1;
% --> (9) Sheguindah Shale and Stony Mountain Formation
%         (from Azmy and Jin, 2018)
Data.DiagenesisFlag(Data.Formation=="Stony Mountain"&Data.LeadAuthor=="Azmy") = 1;
Data.DiagenesisFlag(Data.Formation=="Sheguindah Shale"&Data.LeadAuthor=="Azmy") = 1;

% --> Liu et al., 2023 (for now, because no Durango value)
Data.DiagenesisFlag(Data.PublicationDOI=="https://doi.org/10.1016/j.gr.2022.10.015") = 1;
% --> Latal et al., 2006 (for now, until SCRIPT IS UPDATE TO EXCLUDE
% RESTRICTED ENVIRONMENTS & PhanSST is recompiled)
Data.DiagenesisFlag(Data.PublicationDOI=="http://doi.org/10.1007/s00531-005-0510-3") = 1;
% --> Until PhanSST is recompiled
Data.DiagenesisFlag(Data.Formation == "Xiazhen" & Data.ProxyValue<-6) = 1;
Data.DiagenesisFlag(Data.Stage == "Katian" & Data.ProxyValue == -6.52) = 1;

% (h) Load assumption files
AssumptionFiles = loadassumptionfiles(Corrections);

% (i) Save
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
     load(filename)
%     % If saved as a global file:
%     load("pHLookup.mat", "pHLookup")

%% PART 4: PARSE DATA & CALCULATE Y/YE

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





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRAVEYARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (1) Loop to assign plate ages to the ExpNo file 
% (This step has already been completed and the file has been resaved)

% LOAD FILES
% Plate model age file
PlateAge = readtable(sprintf('%s/Data/SupplementalData/PlateAgeByExpNo.xlsx',OneDrive));
% Exp no file
ExpNo = readtable(sprintf('%s/Data/SupplementalData/ExpNoByStage.xlsx',OneDrive));
% ADD VARS
ExpNo.ManualPlateAge1 = NaN(90,1);
ExpNo.ManualPlateAge2 = NaN(90,1);
ExpNo.ManualPlateAge3 = NaN(90,1);
% RUN LOOP
for n = 1:90
    idx1 = find(PlateAge.ExperNo == ExpNo.ManualExpNo1(n));
    idx2 = find(PlateAge.ExperNo == ExpNo.ManualExpNo2(n));
    idx3 = find(PlateAge.ExperNo == ExpNo.ManualExpNo3(n));
    ExpNo.ManualPlateAge1(n) = PlateAge.PlateModelAge(idx1);
    ExpNo.ManualPlateAge2(n) = mean(PlateAge.PlateModelAge([idx1,idx2]));
    ExpNo.ManualPlateAge3(n) = mean(PlateAge.PlateModelAge([idx1,idx3]));
end
% SAVE
writetable(ExpNo,'ExpNoByStage.xlsx')




