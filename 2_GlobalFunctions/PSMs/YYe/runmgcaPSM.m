function Ye = runmgcaPSM( ...
    ens, ensmeta, exclude, proxydata, age, stagename, pHglobal, Nens)

% (0) Load necessary files
    ensmat = ens.load;
    load("HadCM3Coordinates.mat","Lat","Lon")
    [Lat, Lon] = meshgrid(Lat,Lon); Lon = Lon(:); Lat = Lat(:);
% (1) Preallocate Matrices:
    % (a) Ye matrix
    Ye = NaN(height(proxydata.mg), Nens);
    % (b) Row index vectors
    rowstos = NaN(height(proxydata.mg), 1);
    rowsso = NaN(height(proxydata.mg), 1);

% (2) Determine rows of variables used to estimate local d18sw
    for ii = 1:height(proxydata.mg)
        rowstos(ii,1) = ensmeta.closestLatLon( "tos", ...
            [proxydata.mg.PaleoLat(ii), proxydata.mg.PaleoLon(ii)], ...
            "exclude", exclude);
        rowsso(ii,1) = ensmeta.closestLatLon( "soanom", ...
            [proxydata.mg.PaleoLat(ii), proxydata.mg.PaleoLon(ii)], ...
            "exclude", exclude);
    end
    rowslatlon = mod(rowstos,ensmeta.lengths(1));

% (3) Extract model prior values from sites
    tos = ensmat(rowstos,:);
    so = ensmat(rowsso,:)+35;
    sitelat = Lat(rowslatlon);
    sitelon = Lon(rowslatlon);

% (4) Estimate omega
    omega = omgph(sitelat,sitelon,proxydata.mg.PaleoWaterDepth);
    
% (6) determine taxon and estimate MgCa at each site, for each prior,
%     and with each iteration
  for ii = 1:height(proxydata.mg)
    if contains(string(proxydata.mg.Taxon2(ii)),'ruber')
        taxon = 'ruber';
    elseif contains(string(proxydata.mg.Taxon2(ii)),'bulloides')
        taxon = 'bulloides';
    elseif contains(string(proxydata.mg.Taxon2(ii)),'sacculifer')
        taxon = 'sacculifer';
    elseif contains(string(proxydata.mg.Taxon2(ii)),'pachy') 
        taxon = 'pachy';  
    elseif contains(string(proxydata.mg.Taxon2(ii)),'incompta') 
        taxon = 'incompta';
    else
        taxon = 'all';
    end  
    lnmg = mgcaPSM_forward(age, tos(ii,:), omega(ii), so(ii,:), ...
        pHglobal, proxydata.mg.CleaningMethod(ii), taxon, 1);
    Ye(ii,:) = median(lnmg,2);
  end  
    
end



