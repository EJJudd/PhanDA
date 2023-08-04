function PaleoCoordinates = makecoordinatelookup(directory,saveoutput,savename,deletefiles)

startingdirectory = cd;
cd(directory)
filenames = dir('*.shp');filenames = extractfield(filenames,'name');

if nargin==1
    saveoutput=0;
end

if nargin==3
    deletefiles=0;
end


for i = 1:length(filenames)
    %read in shapefile
    s = struct2table(shaperead(filenames{i}));
    %create fieldname for PaleoCoordinates structure
    fn=strsplit(filenames{i},'_'); fn=strsplit(fn{2},'Ma.shp');
    fn = str2double(fn{1});fn=sprintf('T%3.1fMa',fn(1));
    fn = strrep(fn,'.','_');
    %fill structure
    PaleoCoordinates.(fn)(:,1) = s.Y;
    PaleoCoordinates.(fn)(:,2) = s.X;
end

if saveoutput == 1
    save(savename,'PaleoCoordinates')
end

if deletefiles == 1
    delete *.shx *.prj *.dbf *.shp
end

cd(startingdirectory);
end