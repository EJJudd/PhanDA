%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Script to extract coordinate data from shape files %%%%%%%%%%%%
%%%%%%%%%%%%%% and make lookup table for discrete entries %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Script created: August 2020 (E. Judd)
%Last modified: August 2020 (E. Judd)

% Set directory name
directory = '/Users/emilyjudd/OneDrive - Syracuse University/PhanTASTIC/GPlates/Outputs/DatasetCoordinatesAll';
% Set file save name (if you want to save the lookup table)
savename = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Code/DataAssimilation/4_GlobalFiles/PaleoCoordinates.mat';
PaleoCoordinates = makecoordinatelookup(directory,1,savename);


