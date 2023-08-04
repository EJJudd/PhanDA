function [Ye, localassumptions] = rund18amacroPSM(stagename, proxydata, ...
    SeawaterLookup, Assumptions, Corrections, Preferences, Nens)

% (0) Pretreat data
d18Oswglobal = Assumptions.global.d18Osw.raw;

% (1) Preallocate Matrices:
d18Oswlocal = NaN(height(proxydata.d18a), Nens);
tos = NaN(height(proxydata.d18a), Nens);
so = NaN(height(proxydata.d18a), Nens);
 
% (2) Use lookuptable to find seawater & pH values for each site and prior
for ii = 1:height(proxydata.d18a)
    % d18O Local
    idx = SeawaterLookup.PaleoLat == proxydata.d18a.PaleoLat(ii) & ...
        SeawaterLookup.PaleoLon == proxydata.d18a.PaleoLon(ii) & ...
        contains(SeawaterLookup.Stage, stagename);
    d18Oswlocal(ii,:) = mean(SeawaterLookup.d18Oswlocal(idx,:),1);
    tos(ii,:) = mean(SeawaterLookup.tos(idx,:),1);
    so(ii,:) = mean(SeawaterLookup.so(idx,:),1);
end
d18Oswfinal = d18Oswlocal + d18Oswglobal;
    
% (3) Calculate Ye    
Ye = temp2d18O(tos,Preferences.d18a.eq,Preferences.d18a.min,d18Oswfinal);

% (4) Apply pH correction if selected
d18OpHlocal_ens = zeros(size(d18Oswlocal));
d18OpHlocal_rec = zeros(size(d18Oswlocal));
if Corrections.phcorrection
    for ii = 1:height(proxydata.d18a)
    for jj = 1:Nens %Consider vectorizing later...
    % Calculate using ens value
    d18opH = pHCorrect(repmat(Ye(ii,jj),1000,1), tos(ii,jj), ...
        so(ii,jj), Assumptions.global.pHens(jj), 0);
    d18OpHlocal_ens(ii,jj) = median(d18opH)-Ye(ii,jj);
    % Calculate using reconstructed value
    d18opH = pHCorrect(repmat(Ye(ii,jj),1000,1), tos(ii,jj), ...
        so(ii,jj), Assumptions.global.pHrec, 0);
    d18OpHlocal_rec(ii,jj) = median(d18opH)-Ye(ii,jj);
    end
    end
 end
  
% (5) Finalize outputs
localassumptions.d18Oswlocal = d18Oswlocal;
localassumptions.d18OpHens = d18OpHlocal_ens;
localassumptions.d18OpHrec = d18OpHlocal_rec;

end


