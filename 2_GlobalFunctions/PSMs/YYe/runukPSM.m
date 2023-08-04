function Ye = runukPSM(ens, ensmeta, exclude, proxydata)

% (1) Preallocate Matrices:
    % (a) Row index vectors
    rows_tos = NaN(height(proxydata.uk), 1);
    % (b) Prior values at proxy site matrices
    tos = NaN(height(proxydata.uk), ensmeta.nMembers);

% (2) Determine rows of variables used to estimate local d18sw
    for ii = 1:height(proxydata.uk)
        rows_tos(ii,1) = ensmeta.closestLatLon( "tos", ...
            [proxydata.uk.PaleoLat(ii), proxydata.uk.PaleoLon(ii)], ...
            'exclude', exclude);
    end
    rows_latlon = mod(rows_tos,ensmeta.lengths(1));
    ens_seasonaltos = load(ens.useVariables("tosmm"));

% (3) Extract (seasonal) model prior values from sites
    for ii = 1:height(proxydata.uk)
        tosidx = rows_latlon(ii):ensmeta.lengths(1):ensmeta.lengths(1)*12;
        if strcmpi(string(proxydata.uk.ContinentOcean(ii)),'pa') && proxydata.uk.PaleoLat(ii) > 30
        tos(ii,:) = mean(ens_seasonaltos(tosidx(6:8),:));   
        elseif strcmpi(string(proxydata.uk.ContinentOcean(ii)),'at') && proxydata.uk.PaleoLat(ii) > 30
        tos(ii,:) = mean(ens_seasonaltos(tosidx(8:10),:));   
        elseif strcmpi(string(proxydata.uk.ContinentOcean(ii)),'me')
        tos(ii,:) = mean(ens_seasonaltos(tosidx([1:5 11:12]),:));
        else
        tos(ii,:) = mean(ens_seasonaltos(tosidx,:));
        end
    end

% (4) Estimate uk value and draw Niter values  
    uk = ukPSM_forward(tos(:));
    uk = median(uk,2);
    Ye = reshape(uk,height(proxydata.uk), ensmeta.nMembers);

end
