function Ye_updated = runukPSM_DAresults(ens, ensmeta, proxydata, output, idx)

    tosmm_updated = ensmeta.regrid("tosmm", output.Amean, 'order', ["lon","lat"]);
    proxydata.uk = proxydata.uk(idx,:);
    
    exclude = mean(load(ens.useVariables('tos')),2);
    exclude(~isnan(exclude)) = false;
    exclude(isnan(exclude)) = true;
    exclude = logical(exclude);
    
% % (1) Preallocate Matrices:
    % (a) Row index vectors
    rows_tos = NaN(height(proxydata.uk), 1);
    % (b) Prior values at proxy site matrices
    tos = NaN(height(proxydata.uk), 1);

% (2) Determine rows of variables use
    for ii = 1:height(proxydata.uk)
        rows_tos(ii,1) = ensmeta.closestLatLon( "tos", ...
            [proxydata.uk.PaleoLat(ii), proxydata.uk.PaleoLon(ii)], ...
            'exclude', exclude);
    end
    rows_latlon = mod(rows_tos,ensmeta.lengths(1));
    [sitelon,sitelat] = ind2sub([96,73],rows_latlon);

% (3) Cycle through model prios extract values from sites
    for ii = 1:height(proxydata.uk)
        if strcmpi(string(proxydata.uk.ContinentOcean(ii)),'pa') && proxydata.uk.PaleoLat(ii) > 30
        tos(ii,:) = mean(tosmm_updated(sitelon(ii),sitelat(ii),6:8));   
        elseif strcmpi(string(proxydata.uk.ContinentOcean(ii)),'at') && proxydata.uk.PaleoLat(ii) > 30
        tos(ii,:) = mean(tosmm_updated(sitelon(ii),sitelat(ii),8:10));   
        elseif strcmpi(string(proxydata.uk.ContinentOcean(ii)),'me')
        tos(ii,:) = mean(tosmm_updated(sitelon(ii),sitelat(ii),[1:5 11:12]));
        else
        tos(ii,:) = mean(tosmm_updated(sitelon(ii),sitelat(ii),:));
        end
    end

% (4) Estimate uk value and draw Niter values 
    uk = ukPSM_forward(tos(:));
    Ye_updated = median(uk,2);
    
    
end

