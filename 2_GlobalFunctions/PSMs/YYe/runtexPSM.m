function Ye = runtexPSM(ens, ensmeta, exclude, proxydata)

% (0) Load necessary files
    ensmat = ens.load;
    % lat/lon values (arbitrary, as not used in analog mode)
    lat = 0; lon = 0;
    
% (1) Preallocate Matrices:
    % (a) Ye matrix
    Ye = NaN(height(proxydata.tex), ensmeta.nMembers);
    % (b) Row index vectors
    rows_tos = NaN(height(proxydata.tex), 1);

% (2) Determine rows of variables used to estimate local d18sw
    for ii = 1:height(proxydata.tex)
        rows_tos(ii,1) = min(ensmeta.closestLatLon("tos", ...
            [proxydata.tex.PaleoLat(ii), proxydata.tex.PaleoLon(ii)], ...
            "exclude", exclude)); 
    end

% (3) Extract model prior values from sites
    tos = ensmat(rows_tos,:);

% (4) Estimate tex value and draw Niter values
    load('tex/Output_SpatAg_SST/params_analog',...
        'alpha_samples','beta_samples','tau2_samples');
    load('Data_Input_SpatAg_SST','Data_Input');    
    for ii = 1:height(proxydata.tex)
        for jj = 1:ensmeta.nMembers
            try
            tex = texPSM_forward(lat, lon, tos(ii,jj), 'SST', 'analog', 5, ...
                alpha_samples,beta_samples,tau2_samples,Data_Input);
            catch
            tex = texPSM_forward(lat, lon, tos(ii,jj), 'SST', 'analog', 7.5, ...
                alpha_samples,beta_samples,tau2_samples,Data_Input);
            end
            Ye(ii,jj) = median(tex,2);
        end
    end
end
