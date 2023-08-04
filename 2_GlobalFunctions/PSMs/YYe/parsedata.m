function [proxydata,uniqueproxydata] = parsedata( ...
    Data, Preferences, ageinfo)

% Apply corrections/user defined preferences
 % (1) Parse data by age range or stage
if strcmpi(Preferences.agemethod,'agerange')
    Data = Data(Data.Age>=ageinfo(1) & Data.Age<ageinfo(2),:);
elseif strcmpi(Preferences.agemethod,'stage')
    Data = Data(contains(Data.Stage, ageinfo),:);
else
    warning('Invalid age parsing method provided')
end

% (2) Include or exclude diagenetically altered data
if strcmpi(Preferences.diagenesisflag,'exclude')
    Data(Data.DiagenesisFlag==1,:) = [];
    Data(strcmpi(string(Data.ProxyType),'tex')& ...
        (Data.dRI>0.5 | Data.MI>0.5 | Data.BIT>0.5),:) = [];
    Data(strcmpi(string(Data.ProxyType),'uk')&Data.PaleoLat>70,:) = [];
elseif ~strcmpi(Preferences.diagenesisflag,"exclude") && ....
        ~strcmpi(Preferences.diagenesisflag,"include")
    warning('Invalid diagenesis preference provided')
end

% (3) Remove data with CAI above prespecified value
Data(Data.MaximumCAI>Preferences.maxCAI,:) = [];

% (4) Include ("all") or exclude ("surfaceonly") subsurface taxa
if strcmpi(Preferences.layer,"surfaceonly")
    Data(contains(Data.Taxon3, 'inoceramus', 'IgnoreCase', true),:) = [];
    Data(contains(Data.Ecology, 'thermocline'),:) = [];
    Data(contains(Data.Ecology, 'other') & strcmpi(Data.Taxon1,'pf'),:) = [];
elseif ~strcmpi(Preferences.layer,"surfaceonly") && ~strcmpi(Preferences.layer,"all")
    warning('Invalid depth layer preference provided')
end

% (5) Remove data from non marine environments
for ii = 1:numel(Preferences.environment)
    Data(contains(Data.Environment, Preferences.environment(ii), 'IgnoreCase', true),:) = [];
end

% Cycle through proxies
proxies = unique(Data.ProxyType);

for ii = 1:numel(proxies)
    proxydata.(proxies{ii}) = Data(strcmp(Data.ProxyType,...
        string(proxies(ii))),:);
    % Isolate unique sampling sites
    % if foram data, parse by taxon
    if strcmpi(string(proxies(ii)),'mg') || strcmpi(string(proxies(ii)),'d18cforam')
        uniquetaxa = unique([proxydata.(proxies{ii}).PaleoLat, ...
            proxydata.(proxies{ii}).PaleoLon, ...
            string(proxydata.(proxies{ii}).Taxon1), ...
            string(proxydata.(proxies{ii}).Taxon3)],'rows');
        idx = (1:numel(uniquetaxa(:,4)))';
        idx(~strcmpi(string(uniquetaxa(:,3)),'pf') | ...
                contains(string(uniquetaxa(:,4)),'ruber') | ...
                contains(string(uniquetaxa(:,4)),'bulloides') | ...
                contains(string(uniquetaxa(:,4)),'sacculifer') | ...
                contains(string(uniquetaxa(:,4)),'pachy') | ...
                contains(string(uniquetaxa(:,4)),'incompta')) = [];
        uniquetaxa(idx,4) = 'mixed';
        uniquetaxa = unique(uniquetaxa,'rows');
        uniquecoords = str2double(uniquetaxa(:,1:2)); uniquetaxa(:,1:2) = [];
    
    % if macrofossil data, parse by taxon
    elseif strcmpi(string(proxies(ii)),'d18cmacro')
        uniquetaxa = unique([proxydata.(proxies{ii}).PaleoLat,proxydata.(proxies{ii}).PaleoLon,string(proxydata.(proxies{ii}).Taxon1),string(proxydata.(proxies{ii}).Taxon2)],'rows');
        uniquecoords = str2double(uniquetaxa(:,1:2)); uniquetaxa(:,1:2) = [];
    elseif strcmpi(string(proxies(ii)),'d18a')
        uniquetaxa = unique([proxydata.(proxies{ii}).PaleoLat,proxydata.(proxies{ii}).PaleoLon,string(proxydata.(proxies{ii}).Taxon1),string(proxydata.(proxies{ii}).Taxon2)],'rows');
        uniquecoords = str2double(uniquetaxa(:,1:2)); uniquetaxa(:,1:2) = [];       
    
    % if d18Op data, parse by methodology
    elseif strcmpi(string(proxies(ii)),'d18p')
        uniquedata = unique([proxydata.(proxies{ii}).PaleoLat,proxydata.(proxies{ii}).PaleoLon,string(proxydata.(proxies{ii}).AnalyticalTechnique)],'rows');
        uniquecoords = str2double(uniquedata(:,1:2)); uniquedata(:,1:2) = [];
    
    % for all other data, parse just by coordinates
    else 
        uniquecoords = unique([proxydata.(proxies{ii}).PaleoLat,proxydata.(proxies{ii}).PaleoLon],'rows');
    end
    
    % Cycle through sampling sites and average data
    pat = ["ruber","bulloides","sacculifer","pachy","incompta"];
    for jj = 1:numel(uniquecoords)/2
    % if foram data, parse by taxon
    if strcmpi(string(proxies(ii)),'mg') || strcmpi(string(proxies(ii)),'d18cforam')
        if strcmpi(uniquetaxa(jj,1),'pf') && ~contains(uniquetaxa(jj,2),pat)
        idx = find(proxydata.(proxies{ii}).PaleoLat ==uniquecoords(jj,1) & ...
           proxydata.(proxies{ii}).PaleoLon ==uniquecoords(jj,2) & ...
            strcmpi(string(proxydata.(proxies{ii}).Taxon1),uniquetaxa(jj,1)) & ...
            ~contains(string(proxydata.(proxies{ii}).Taxon3),pat));
        else
        idx = find(proxydata.(proxies{ii}).PaleoLat ==uniquecoords(jj,1) & ...
           proxydata.(proxies{ii}).PaleoLon ==uniquecoords(jj,2) & ...
            strcmpi(string(proxydata.(proxies{ii}).Taxon1),uniquetaxa(jj,1)) & ...
            contains(string(proxydata.(proxies{ii}).Taxon3),uniquetaxa(jj,2)));
        end
    
    % if macrofossil data, parse by taxon
    elseif strcmpi(string(proxies(ii)),'d18cmacro') || strcmpi(string(proxies(ii)),'d18a')
       idx = find(proxydata.(proxies{ii}).PaleoLat ==uniquecoords(jj,1) & ...
        proxydata.(proxies{ii}).PaleoLon ==uniquecoords(jj,2) & ...
        strcmpi(string(proxydata.(proxies{ii}).Taxon1),uniquetaxa(jj,1)) & ...
        contains(string(proxydata.(proxies{ii}).Taxon2),uniquetaxa(jj,2)));
    
    % if d18p data, parse by methodology
    elseif strcmpi(string(proxies(ii)),'d18p')
       idx = find(proxydata.(proxies{ii}).PaleoLat ==uniquecoords(jj,1) & ...
        proxydata.(proxies{ii}).PaleoLon ==uniquecoords(jj,2) & ...
        strcmpi(string(proxydata.(proxies{ii}).AnalyticalTechnique),uniquedata(jj,1)));
    
    % if any other proxy, just find unique coords
    else 
       idx = find(proxydata.(proxies{ii}).PaleoLat ==uniquecoords(jj,1) & ...
         proxydata.(proxies{ii}).PaleoLon ==uniquecoords(jj,2));
    end
    
    uniqueproxydata.(proxies{ii}).SiteName(jj,1) = join(unique(proxydata.(proxies{ii}).SiteName(idx)),', ');
    uniqueproxydata.(proxies{ii}).Age(jj,1) = mean(proxydata.(proxies{ii}).Age(idx));
    uniqueproxydata.(proxies{ii}).ModLat(jj,1) = mean(proxydata.(proxies{ii}).ModLat(idx));
    uniqueproxydata.(proxies{ii}).ModLon(jj,1) = mean(proxydata.(proxies{ii}).ModLon(idx));
    uniqueproxydata.(proxies{ii}).PaleoLat(jj,1) = unique(proxydata.(proxies{ii}).PaleoLat(idx));
    uniqueproxydata.(proxies{ii}).PaleoLon(jj,1) = unique(proxydata.(proxies{ii}).PaleoLon(idx));
    uniqueproxydata.(proxies{ii}).ProxyValue(jj,1) = median(proxydata.(proxies{ii}).ProxyValue(idx));
    uniqueproxydata.(proxies{ii}).ProxyVar(jj,1) = var(proxydata.(proxies{ii}).ProxyValue(idx));
    uniqueproxydata.(proxies{ii}).N(jj,1) = numel(idx);
    uniqueproxydata.(proxies{ii}).AnalyticalMethod(jj,1) = join(unique(proxydata.(proxies{ii}).AnalyticalTechnique(idx)),'; ');
    uniqueproxydata.(proxies{ii}).ContinentOcean(jj,1) = join(unique(proxydata.(proxies{ii}).ContinentOcean(idx)),', ');
    if strcmpi(string(proxies(ii)),'mg') || contains(string(proxies(ii)),'d18c')
    uniqueproxydata.(proxies{ii}).Taxon1(jj,1) = uniquetaxa(jj,1);
    uniqueproxydata.(proxies{ii}).Taxon2(jj,1) = uniquetaxa(jj,2);
    else
    uniqueproxydata.(proxies{ii}).Taxon1(jj,1) = join(unique(proxydata.(proxies{ii}).Taxon1(idx)),', ');
    uniqueproxydata.(proxies{ii}).Taxon2(jj,1) = join(unique(proxydata.(proxies{ii}).Taxon2(idx)),', ');    
    end
    if strcmpi(string(proxies(ii)),'mg')
    uniqueproxydata.(proxies{ii}).lnmg(jj,1) = log(uniqueproxydata.(proxies{ii}).ProxyValue(jj,1));
    uniqueproxydata.(proxies{ii}).CleaningMethod(jj,1) = mean(proxydata.(proxies{ii}).CleaningMethod(idx));
    uniqueproxydata.(proxies{ii}).PaleoWaterDepth(jj,1) = nanmean(proxydata.(proxies{ii}).PaleoWaterDepth(idx));
    else
    uniqueproxydata.(proxies{ii}).lnmg(jj,1) = NaN;
    uniqueproxydata.(proxies{ii}).CleaningMethod(jj,1) = NaN;
    uniqueproxydata.(proxies{ii}).PaleoWaterDepth(jj,1) = NaN;
    end
    authorlist = join([string(proxydata.(proxies{ii}).LeadAuthor(idx)), num2str(proxydata.(proxies{ii}).Year(idx))],', ');
    uniqueproxydata.(proxies{ii}).Ref(jj,1) = join(unique(authorlist),'; ');
    end
    uniqueproxydata.(proxies{ii}) = struct2table(uniqueproxydata.(proxies{ii}));
end

end
