function [Ye, localassumptions] = rund18cforamPSM(stagename, proxydata, ...
    SeawaterLookup, Assumptions, Corrections, Nens)
    
% (0) Pretreat data - isolate just foram data
d18Oswglobal = Assumptions.global.d18Osw.raw;
proxydata.d18c = proxydata.d18cforam;
    
% (1) Preallocate Matrices:
Ye = NaN(height(proxydata.d18c), Nens);
d18Oswlocal = NaN(height(proxydata.d18c), Nens);
tos = NaN(height(proxydata.d18c), Nens);
so = NaN(height(proxydata.d18c), Nens);

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

% (3) determine taxon and estimate d18foram at each site, for each prior,
%     and with each iteration
  for ii = 1:height(proxydata.d18c)
    if contains(string(proxydata.d18c.Taxon2(ii)),'ruber')
        taxon = 'ruber';
    elseif contains(string(proxydata.d18c.Taxon2(ii)),'bulloides')
        taxon = 'bulloides';
    elseif contains(string(proxydata.d18c.Taxon2(ii)),'sacculifer')
        taxon = 'sacculifer';
    elseif contains(string(proxydata.d18c.Taxon2(ii)),'pachy') 
        taxon = 'pachy';  
    elseif contains(string(proxydata.d18c.Taxon2(ii)),'incompta') 
        taxon = 'incompta';
    else
        taxon = 'all';
    end  
    d18foramEst = d18foram_forward(tos(ii,:)', d18Oswfinal(ii,:)', taxon);
    Ye(ii,:) = median(d18foramEst,2);
  end

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


