function mapfilenames = selectmap(timeplot, mapdirectory)
%Timeplot: vector of ages you want to plot
%Mapdirectory: where your map csvs are stores
startingdirectory = cd;
cd(mapdirectory)
filenames = dir('*.csv');
filenames = extractfield(filenames,'name')';

for ii = 1:length(filenames)
    fn=strsplit(filenames{ii},'_');fn=strsplit(fn{2},'Ma');
    MapTime(ii,1) = str2num(fn{1});
end

mapfilenames = cell(numel(timeplot),1);
for ii = 1:numel(timeplot)
    mapselect = filenames(abs(MapTime-timeplot(ii))==min(abs(MapTime-timeplot(ii))));
    mapfilenames{ii} = sprintf('%s/%s', mapdirectory, string(mapselect{1}));
end

cd(startingdirectory)

end