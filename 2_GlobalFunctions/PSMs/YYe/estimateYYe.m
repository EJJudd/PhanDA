function [Y, Ye, Assumptions] = estimateYYe(proxydata, ens, ensmeta, ...
    stageno, stagename, SeawaterLookup, pHLookup, AssumptionFiles, ...
    Corrections, Preferences, GTS, verbose)

% (0) Pretreat/load values
    % (a) Proxy types
    fn = fieldnames(proxydata);
    % (b) Number of ensemble members
    Nens = ensmeta.nMembers;
    % (c) Location of data to exclude in closestLatLon function
    exclude = mean(load(ens.useVariables('tos')),2);
    exclude(~isnan(exclude)) = false;
    exclude(isnan(exclude)) = true;
    exclude = logical(exclude);

% (1) Generate secular assumptions
     load("PhanerozoicpHv6.mat","PhanerozoicpH")
     % (a) Load seawater curves 
     Assumptions.global.pHens = pHLookup.pHGlobal(stageno,:);
     Assumptions.global.pHrec = median(PhanerozoicpH(stageno,:));
    % (b) Estimate Global d18Osw
    Assumptions.global.d18Osw.raw = AssumptionFiles.Globald18O(stageno);
    % (c) Estimate snowball earth correction
    Assumptions.global.d18Osw.Snowball = AssumptionFiles.Snowball(stageno);

% (2) Loop through proxies and assemble Y, Ye values
for ii = 1:numel(fn)

    % Assign Y value 
    % (if Mg/Ca, this will get replace below with the ln value)
    Y.(fn{ii}) = proxydata.(fn{ii}).ProxyValue;
    
    % Initialize printprogress messages
    message1 = sprintf('Beginning proxy %.0f of %.0f (%s)',ii,numel(fn),fn{ii});
    message2 = sprintf('Completed proxy %.0f of %.0f (%s)',ii,numel(fn),fn{ii});
    
    % Assign Ye values
    if strcmpi(fn(ii),'d18a')
        printprogress(message1, verbose)
        [Ye.d18a, Assumptions.local.d18a] = rund18amacroPSM(stagename, ...
            proxydata, SeawaterLookup, Assumptions, Corrections, ...
            Preferences, Nens);
        printprogress(message2, verbose)
    
    elseif strcmpi(fn(ii),'d18cmacro')
        printprogress(message1, verbose)
        [Ye.d18cmacro, Assumptions.local.d18cmacro] = rund18cmacroPSM(...
            stagename, proxydata, SeawaterLookup, Assumptions, ...
            Corrections, Preferences, Nens);
        printprogress(message2, verbose)
    
    elseif strcmpi(fn(ii),'d18cforam')
        printprogress(message1, verbose)
        [Ye.d18cforam, Assumptions.local.d18cforam] = rund18cforamPSM(...
            stagename, proxydata, SeawaterLookup, Assumptions, ...
            Corrections, Nens);      
        printprogress(message2, verbose)
    
    elseif strcmpi(fn(ii),'d18p')
        printprogress(message1, verbose)
        [Ye.d18p, Assumptions.local.d18p] = rund18pPSM(stagename, ...
            proxydata, SeawaterLookup, Assumptions, Preferences, Nens);    
        printprogress(message2, verbose)
        
    elseif contains(fn(ii),'mg')
        printprogress(message1, verbose)
        Y.mg = proxydata.mg.lnmg;
        Ye.mg = runmgcaPSM(ens, ensmeta, exclude, proxydata, ...
            GTS.Average(stageno), stagename, Assumptions.global.pH, Nens);
        printprogress(message2, verbose)
    
    elseif strcmpi(fn(ii),'tex')
        printprogress(message1, verbose)
        Ye.tex = runtexPSM(ens, ensmeta, exclude, proxydata);
        printprogress(message2, verbose)
    
    elseif strcmpi(fn(ii),'uk')
        printprogress(message1, verbose)
        Ye.uk = runukPSM(ens, ensmeta, exclude, proxydata);
        printprogress(message2, verbose)
    
    end
    
end

if ~exist('Y')
    Y = [];
    Ye = [];
end

end