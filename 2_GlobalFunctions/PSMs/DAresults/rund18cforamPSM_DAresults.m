function Ye_updated = ...
    rund18cforamPSM_DAresults(ens, ensmeta, ensmat, proxydata, output, d18swglobal, ...
    pHglobal, AssumptionFiles, idx)
    
    load("HadCM3Coordinates.mat","Lat","Lon")
    [Lat, ~] = meshgrid(Lat,Lon);

    tos_updated = ensmeta.regrid("tos", output.Amean, 'order', ["lon","lat"]);
    so_updated = 35 + ensmeta.regrid("soanom", output.Amean, 'order', ["lon","lat"]);
    tosanom_updated = ensmeta.regrid("tosanom", output.Amean, 'order', ["lon","lat"]);
    soanom_updated = ensmeta.regrid("soanom", output.Amean, 'order', ["lon","lat"]);
    pranom_updated = ensmeta.regrid("prlnanom", output.Amean, 'order', ["lon","lat"]);
    dist_updated = ensmeta.regrid("dist", ensmat, 'order', ["lon","lat"]);
    proxydata.d18cforam = proxydata.d18cforam(idx,:);

    exclude = mean(load(ens.useVariables('tos')),2);
    exclude(~isnan(exclude)) = false;
    exclude(isnan(exclude)) = true;
    exclude = logical(exclude);
    
% (1) Preallocate Matrices:
    Ye_updated = NaN(height(proxydata.d18cforam), 1);
    rows_tos = NaN(height(proxydata.d18cforam), 1);

% (2) Determine rows of variables use
    for ii = 1:height(proxydata.d18cforam)
        rows_tos(ii,1) = ensmeta.closestLatLon( "tos", ...
            [proxydata.d18cforam.PaleoLat(ii), proxydata.d18cforam.PaleoLon(ii)], ...
            'exclude', exclude);
    end
    rows_latlon = mod(rows_tos,ensmeta.lengths(1));
    tos = tos_updated(rows_latlon);
    so = so_updated(rows_latlon);
    x.sst_anomaly = tosanom_updated(rows_latlon);
    x.sss_anomaly = soanom_updated(rows_latlon);
    x.precipln_anomaly = pranom_updated(rows_latlon);
    x.dist = dist_updated(rows_latlon);
    x.lat = Lat(rows_latlon);
    
% (3) Estimate local seawater
    d18swlocal = predictd18Osw(x, AssumptionFiles, true);
    d18swfinal = d18swlocal + d18swglobal;
    
% (6) determine taxon and estimate d18foram at each site, for each prior,
%     and with each iteration
  for ii = 1:height(proxydata.d18cforam)
    if contains(string(proxydata.d18cforam.Taxon2(ii)),'ruber')
        taxon = 'ruber';
    elseif contains(string(proxydata.d18cforam.Taxon2(ii)),'bulloides')
        taxon = 'bulloides';
    elseif contains(string(proxydata.d18cforam.Taxon2(ii)),'sacculifer')
        taxon = 'sacculifer';
    elseif contains(string(proxydata.d18cforam.Taxon2(ii)),'pachy') 
        taxon = 'pachy';  
    elseif contains(string(proxydata.d18cforam.Taxon2(ii)),'incompta') 
        taxon = 'incompta';
    else
        taxon = 'all';
    end  
    d18foram = d18foram_forward(tos(ii), d18swfinal(ii), taxon);
    d18foramph = pHCorrect(d18foram,tos(ii),so(ii),pHglobal);
    Ye_updated(ii,1) = median(d18foramph);
  end


end


