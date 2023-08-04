function Ye_updated = runtexPSM_DAresults(ens, ensmeta, proxydata, output, idx)

    tos_updated = ensmeta.regrid("tos", output.Amean, 'order', ["lon","lat"]);
    proxydata.tex = proxydata.tex(idx,:);

    exclude = mean(load(ens.useVariables('tos')),2);
    exclude(~isnan(exclude)) = false;
    exclude(isnan(exclude)) = true;
    exclude = logical(exclude);
    

% (1) Preallocate Matrices:
    % (a) Ye matrix
    Ye_updated = NaN(height(proxydata.tex), 1);
    % (b) Row index vectors
    rows_tos = NaN(height(proxydata.tex), 1);
    % (c) Prior values at proxy site matrices
    tos = NaN(height(proxydata.tex),1);

% (2) Determine rows of variables use
    for ii = 1:height(proxydata.tex)
        rows_tos(ii,1) = ensmeta.closestLatLon( "tos", ...
            [proxydata.tex.PaleoLat(ii), proxydata.tex.PaleoLon(ii)], ...
            'exclude', exclude);
    end
    rows_latlon = mod(rows_tos,ensmeta.lengths(1));
    tos(:,1) = tos_updated(rows_latlon);
   
% (4) Estimate tex value and draw Niter values 
    load('tex/Output_SpatAg_SST/params_analog',...
        'alpha_samples','beta_samples','tau2_samples');
    load('Data_Input_SpatAg_SST','Data_Input');    
    for ii = 1:height(proxydata.tex)
        try
        texEst = texPSM_forward(0, 0, tos(ii), ...
            'SST', 'analog', 7.5, alpha_samples, beta_samples, ...
            tau2_samples, Data_Input);
        catch
        texEst = texPSM_forward(0,0, tos(ii), ...
            'SST', 'analog', 15, alpha_samples, beta_samples, ...
            tau2_samples, Data_Input);
        end
        Ye_updated(ii,1) = mean(texEst);
    end


end
