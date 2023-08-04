function [Ye, localassumptions] = rund18pPSM(stagename, proxydata, ...
    SeawaterLookup, Assumptions, Preferences, Nens)

% (0) Pretreat data
d18Oswglobal = Assumptions.global.d18Osw.raw;

% (1) Preallocate Matrices:
d18Oswlocal = NaN(height(proxydata.d18p), Nens);
tos = NaN(height(proxydata.d18p), Nens);

% (2) Use lookuptable to find seawater value for each site and prior
for ii = 1:height(proxydata.d18p)
    idx = find(SeawaterLookup.PaleoLat == proxydata.d18p.PaleoLat(ii) & ...
        SeawaterLookup.PaleoLon == proxydata.d18p.PaleoLon(ii) & ...
        contains(SeawaterLookup.Stage, stagename));
    d18Oswlocal(ii,:) = mean(SeawaterLookup.d18Oswlocal(idx,:),1);
    tos(ii,:) = mean(SeawaterLookup.tos(idx,:),1);    
end
d18Oswfinal = d18Oswlocal + d18Oswglobal;
    
% (3) Calculate Ye    
Ye = temp2d18O(tos,Preferences.d18p.eq,Preferences.d18p.min,d18Oswfinal);
    
% (6) Finalize outputs
localassumptions.d18Oswlocal = d18Oswlocal;

end



