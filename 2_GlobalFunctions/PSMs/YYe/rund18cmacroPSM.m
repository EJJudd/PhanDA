function [Ye, localassumptions] = rund18cmacroPSM(stagename, proxydata, ...
    SeawaterLookup, Assumptions, Corrections, Preferences, Nens)

% (0) Pretreat data - isolate just macrofossils
d18Oswglobal = Assumptions.global.d18Osw.raw;
proxydata.d18c = proxydata.d18cmacro;

% (1) Preallocate Matrices:
d18Oswlocal = NaN(height(proxydata.d18c), Nens);
tos = NaN(height(proxydata.d18c), Nens);
so = NaN(height(proxydata.d18c), Nens);
Ye = NaN(height(proxydata.d18c), Nens);

% (2) Use lookuptable to find seawater value for each site and prior
for ii = 1:height(proxydata.d18c)
    % d18O Local
    idx = SeawaterLookup.PaleoLat == proxydata.d18c.PaleoLat(ii) & ...
        SeawaterLookup.PaleoLon == proxydata.d18c.PaleoLon(ii) & ...
        contains(SeawaterLookup.Stage, stagename);
    d18Oswlocal(ii,:) = mean(SeawaterLookup.d18Oswlocal(idx,:),1);
    tos(ii,:) = mean(SeawaterLookup.tos(idx,:),1);
    so(ii,:) = mean(SeawaterLookup.so(idx,:),1);
end
d18Oswfinal = d18Oswlocal + d18Oswglobal;

% (3) Calculate Ye   
nobelidx = find(~strcmpi(proxydata.d18c.Taxon2,'ce'));
belidx =  find(strcmpi(proxydata.d18c.Taxon2,'ce'));
Ye(nobelidx,:) = temp2d18O(tos(nobelidx,:), Preferences.d18c.eq.nobel,...
    Preferences.d18c.min, d18Oswfinal(nobelidx,:));
Ye(belidx,:) = temp2d18O(tos(belidx,:), Preferences.d18c.eq.bel,...
    Preferences.d18c.min, d18Oswfinal(belidx,:));

% (4) Apply pH correction if selected
d18OpHlocal_ens = zeros(size(d18Oswlocal));
d18OpHlocal_rec = zeros(size(d18Oswlocal));
if Corrections.phcorrection
    for ii = 1:height(proxydata.d18c)
    for jj = 1:Nens %Consider vectorizing later...
    % Calculate using ens value
    d18opH = pHCorrect(repmat(Ye(ii,jj),1000,1), tos(ii,jj), ...
        so(ii,jj), Assumptions.global.pHens(jj), 0);
    d18OpHlocal_ens(ii,jj) = median(d18opH)-Ye(ii,jj);
    d18opH = pHCorrect(repmat(Ye(ii,jj),1000,1), tos(ii,jj), ...
        so(ii,jj), Assumptions.global.pHrec, 0);
    d18OpHlocal_rec(ii,jj) = median(d18opH)-Ye(ii,jj);
    end
    end
end

    
% (5) Correct for belemenites (does nothing if correction is set to 0)
Ye(strcmpi(proxydata.d18c.Taxon2,"ce"),:) = ...
    Ye(strcmpi(proxydata.d18c.Taxon2,"ce"),:) + Corrections.belemnite;

     
% (6) Finalize outputs
localassumptions.d18Oswlocal = d18Oswlocal;
localassumptions.d18OpHens = d18OpHlocal_ens;
localassumptions.d18OpHrec = d18OpHlocal_rec;

end



