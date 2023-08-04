function [tos,omega,so,pH,taxon] = runmgcaPSM_DAresults(ens, ensmeta, ...
    proxydata, output, idx)

    tos_updated = ensmeta.regrid("tos", output.Amean, 'order', ["lon","lat"]);
    so_updated = 35 + ensmeta.regrid("soanom", output.Amean, 'order', ["lon","lat"]);
    proxydata.mg = proxydata.mg(idx,:);

    exclude = mean(load(ens.useVariables('tos')),2);
    exclude(~isnan(exclude)) = false;
    exclude(isnan(exclude)) = true;
    exclude = logical(exclude);
        
% (1) Preallocate Matrices:
    rows_tos = NaN(height(proxydata.mg), 1);

% (2) Determine rows of variables use
    for ii = 1:height(proxydata.mg)
        rows_tos(ii,1) = ensmeta.closestLatLon( "tos", ...
            [proxydata.mg.PaleoLat(ii), proxydata.mg.PaleoLon(ii)], ...
            'exclude', exclude);
    end
    rows_latlon = mod(rows_tos,ensmeta.lengths(1));
    tos = tos_updated(rows_latlon);
    so = so_updated(rows_latlon);

% (4) Estimate omega
    [omega, pH] = omgph(proxydata.mg.PaleoLat,proxydata.mg.PaleoLon,proxydata.mg.PaleoWaterDepth);
    
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
%     lnmg = mgcaPSM_forward(proxydata.mg.Age, tos(ii), omega(ii), so(ii), ...
%         pH, proxydata.mg.CleaningMethod(ii), taxon, 1);
%     Ye_updated(ii,:) = median(lnmg,2);
  end  
    
end



