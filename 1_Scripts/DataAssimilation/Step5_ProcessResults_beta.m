%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSIMILATION STEP 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process Results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 08/22 (E. Judd)
% Last updated: 08/22 (E. Judd)

% Notes: 

%   -->  This script is compatable with the most recent Dash release
%        (https://github.com/JonKing93/DASH/releases/tag/v4.0.0-beta-2)



% Define dassimilation directories
% (a) Home Directory (root for all other directories)
OneDrive = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC';
iCloud = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs';
% (b) Load the rest of the directories off of the "DirectoryInfo" file stored
% in the InputWorkspaces folder
% If working with assimilation input workspaces made the same day:
% load( sprintf( '%s/AssimilationOutputs/PhanerozoicDA_%s/InputWorkspaces/DirectoryInfo.mat', ...
%                      HomeDir, strrep(date,'-','') ) );
% If working off of older input workspace, manually enter the date
% (Convention: ddMmmyyyy)
filedate = '09Mar2023';
AssDir = [OneDrive,'/AssimilationOutputs/PhanerozoicDA_',filedate];
InputDir = [AssDir,'/InputWorkspaces'];
OutputDir = [AssDir,'/OutputWorkspaces'];
FigDir = [AssDir,'/OutputFigs'];
OutputFilename = 'Output_Rpercentile_SBtrue_CO2true.mat';

% Define additional assimilation & figure parameters
DArange = 1:91;
overwrite = false;

% FIGURE SET 1: BATCH PRINT SUMMARY OUTPUT PLOTS BY STAGE
printsummaryoutputplots(OneDrive, InputDir, OutputDir, FigDir, OutputFilename, DArange, overwrite)

% FIGURE SET 2: PHANEROZOIC GMST CURVES
agerange = [0,GTS.LowerBoundary(DArange(end))];
SaveName = 'PhanGMST_Theory.png';
opts = "line";
printsummarygmst(OutputFilename, OutputDir, FigDir, SaveName, ...
    DArange, agerange, opts)

% FIGURE SET 2: PHANEROZOIC GMST CURVES
agerange = [0,GTS.LowerBoundary(23)];
SaveName = 'Cenoz.png';
opts = "line";
printsummarygmst(OutputFilename, OutputDir, FigDir, SaveName, ...
    DArange, agerange, opts)