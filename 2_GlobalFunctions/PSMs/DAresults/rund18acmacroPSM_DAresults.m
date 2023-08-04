function Ye_updated = ...
    rund18acmacroPSM_DAresults(ens, ensmeta, ensmat, proxydata, output, d18swglobal, ...
    pHglobal, AssumptionFiles, idx)

    load("HadCM3Coordinates.mat","Lat","Lon")
    [Lat, ~] = meshgrid(Lat,Lon);

    tos_updated = ensmeta.regrid("tos", output.Amean, 'order', ["lon","lat"]);
    so_updated = 35 + ensmeta.regrid("soanom", output.Amean, 'order', ["lon","lat"]);
    tosanom_updated = ensmeta.regrid("tosanom", output.Amean, 'order', ["lon","lat"]);
    soanom_updated = ensmeta.regrid("soanom", output.Amean, 'order', ["lon","lat"]);
    pranom_updated = ensmeta.regrid("prlnanom", output.Amean, 'order', ["lon","lat"]);
    dist_updated = ensmeta.regrid("dist", ensmat, 'order', ["lon","lat"]);
    proxydata.d18acmacro = proxydata.d18acmacro(idx,:);

    exclude = mean(load(ens.useVariables('tos')),2);
    exclude(~isnan(exclude)) = false;
    exclude(isnan(exclude)) = true;
    exclude = logical(exclude);

% (1) Preallocate Matrices:
    % (a) Ye matrix
    Ye_updated = NaN(height(proxydata.d18acmacro), 1);
    % (b) Row index vectors
    rows_tos = NaN(height(proxydata.d18acmacro), 1);

% (2) Determine rows of variables used to estimate local d18sw
    for ii = 1:height(proxydata.d18acmacro)
        rows_tos(ii,1) = ensmeta.closestLatLon( "tos", ...
            [proxydata.d18acmacro.PaleoLat(ii), proxydata.d18acmacro.PaleoLon(ii)], ...
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
    for ii = 1:height(proxydata.d18acmacro)
        if proxydata.d18acmacro.ProxyType == "d18cmacro"
            d18macro = temp2d18O(tos(ii),"Kim1997","calcite",d18swfinal);
        elseif proxydata.d18acmacro.ProxyType == "d18a"
            d18macro = temp2d18O(tos(ii),"Grossman1986","aragonite",d18swfinal);
        end
        d18macroph = pHCorrect(repmat(d18macro,1000,1),tos(ii),so(ii),pHglobal);
        Ye_updated(ii,1) = median(d18macroph);
    end
    
end



