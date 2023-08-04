%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Script to extract coordinate data from compilation file %%%%%%%%%%
%%%%%%%%%%%%%%% and make shape file to use with GPlates %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Script created: August 2020 (E. Judd)
%Last modified: August 2020 (E. Judd)


% Extract lat/lon data and create list of unique values
coords(:,1) = Data.ModLat; coords(:,2) = Data.ModLon;
coords = unique(coords,'rows');coords(isnan(coords(:,1)),:) = [];
% Make shape file
filename = '/Users/emilyjudd/OneDrive - Syracuse University/PhanTASTIC/GPlates/Data/DatasetCoordinatesAll.shp';
makeshapefile(coords,filename);