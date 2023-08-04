function pHLookup = makepHlookuptable(DArange, Nens, ExpNo, EnsMat, EnsMeta, ...
    GTS, Lat, Lon, filename)

% (1) Preload necessary information
[Lat,Lon] = meshgrid(Lat,Lon);
Lat = Lat(:); Lon = Lon(:);


dims = string(EnsMeta.ensembleDimensions{1});
expno = EnsMeta.ensemble{1,1}{1, dims == "expno"};
load("ExpCO2.mat","ExpCO2")


% (c) make seperate tables for each ensemble scenario, find unique coors/
pHLookup.pHGlobal = NaN(numel(DArange),Nens);
pHLookup.ExpCO2 = NaN(numel(DArange),Nens);
pHLookup.GMSST = NaN(numel(DArange),Nens);
pHLookup.EnsembleIndex = NaN(numel(DArange),Nens);

% (4) Cycle through the stages to estimate local d18Osw
for ii = DArange
    % (a) index the stage, the sites of that stage, and ensemble members
    sidx = find(ExpNo.Stage == GTS.Stage(ii));
    expidx = [ExpNo.ExpNo1(sidx), ...
              ExpNo.ExpNo2(sidx)];
    ensidx = any(expno == expidx,2);
    pHLookup.EnsembleIndex(ii,:) = find(ensidx);
    ensusemat = EnsMat(:,ensidx);
    ensusemeta = EnsMeta.removeMembers(~ensidx);
    
    tos = ensusemeta.regrid("tos",ensusemat,'order',["lat";"lon"]);
    pHLookup.GMSST(ii,:) = latweightgmst(tos);
            
    exprun = ensusemeta.ensemble{1,1}{1, dims == "exprun"};
    expname = unique(exprun);
    for jj = 1:numel(expname)
        eidx = exprun == expname(jj);
        pHLookup.ExpCO2(ii,eidx) = mean(...
            ExpCO2.(expname{jj})(any(ExpCO2.ExpNo == expidx,2)));
    end
    fprintf('Completed entry for Stage %d/%d\n', ii, height(GTS))

end

% (5) Estimate local pH seawater using CO2SYS
% CO2
co2 = pHLookup.ExpCO2(:);
% Use modern distribution of alk.
alk = normrnd(2300,100,length(co2),1);
% T from pruor
T = pHLookup.GMSST(:);
% Salinity 
S = normrnd(34,2,length(co2),1);
% Define constants for carbonate system calculations
SCALE  = 1; % Total pH scale
K1K2   = 4; % K1 K2 from Mehrbach (1973), refit by Dickson & Millero (1987)
KSO4   = 1; % KSO4 from Dickson (1990) 
KF     = 2; % KHF from Perez & Fraga (1987)
BOR    = 2; % Borate-to-salinity ratio from Lee et al (2010)
% CO2SYS call. 
% Parameter 1 is alkalinity (1). Parameter 2 is CO2 (4). S and T come next,
% then Tout (and all outher 'out' values) is NaN. Surface pressure is 0.
% Silicate, phosphate, ammonium, H2S just set to arbitrary low values. 
A = CO2SYS(alk,co2,1,4,S,T,nan,0,nan,1,1,0,0,SCALE,K1K2,KSO4,KF,BOR);
pH = A(:,3);
pHLookup.pHGlobal = reshape(pH,size(pHLookup.ExpCO2));

% (7) Save file
if exist('filename')
    save(filename,'pHLookup');
end

end
    
    