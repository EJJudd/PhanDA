function Ye_updated = ...
    rund18pPSM_DAresults(ens, ensmeta, ensmat, proxydata, output, d18swglobal, ...
    AssumptionFiles, idx)

    load("HadCM3Coordinates.mat","Lat","Lon")
    [Lat, ~] = meshgrid(Lat,Lon);

    tos_updated = ensmeta.regrid("tos", output.Amean, 'order', ["lon","lat"]);
    so_updated = 35 + ensmeta.regrid("soanom", output.Amean, 'order', ["lon","lat"]);
    tosanom_updated = ensmeta.regrid("tosanom", output.Amean, 'order', ["lon","lat"]);
    soanom_updated = ensmeta.regrid("soanom", output.Amean, 'order', ["lon","lat"]);
    pranom_updated = ensmeta.regrid("prlnanom", output.Amean, 'order', ["lon","lat"]);
    dist_updated = ensmeta.regrid("dist", ensmat, 'order', ["lon","lat"]);
    proxydata.d18p = proxydata.d18p(idx,:);

    exclude = mean(load(ens.useVariables('tos')),2);
    exclude(~isnan(exclude)) = false;
    exclude(isnan(exclude)) = true;
    exclude = logical(exclude);

% (1) Preallocate Matrices:
    % (a) Ye matrix
    Ye_updated = NaN(height(proxydata.d18p), 1);
    % (b) Row index vectors
    rows_tos = NaN(height(proxydata.d18p), 1);

% (2) Determine rows of variables used to estimate local d18sw
    for ii = 1:height(proxydata.d18p)
        rows_tos(ii,1) = ensmeta.closestLatLon( "tos", ...
            [proxydata.d18p.PaleoLat(ii), proxydata.d18p.PaleoLon(ii)], ...
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
    d18swlocal = predictd18Osw(x, AssumptionFiles, false);
    d18swfinal = d18swlocal + d18swglobal;

% (4) Forward model d18O
    Ye_updated = temp2d18O(tos,"Lecuyer2013","phosphate",d18swfinal);
    
end



