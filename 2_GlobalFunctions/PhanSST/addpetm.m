function Data = addpetm(Data, Preferences)

load("PETM.mat","PETM")

if Preferences.PETM == "deepMIP"
    PETM.StartDeepMIP(isnan(PETM.StartDeepMIP)) = ...
        PETM.Start(isnan(PETM.StartDeepMIP));
    PETM.EndDeepMIP(isnan(PETM.EndDeepMIP)) = ...
        PETM.EndPeak(isnan(PETM.EndDeepMIP));
end

for ii = 1:height(PETM)
    if Preferences.PETM == "peak"
        depth(1) = PETM.Start(ii);
        depth(2) = PETM.EndPeak(ii);
    elseif Preferences.PETM == "body"
        depth(1) = PETM.Start(ii);
        depth(2) = PETM.EndPeak(ii); 
    elseif Preferences.PETM == "deepMIP"
        depth(1) = PETM.StartDeepMIP(ii);
        depth(2) = PETM.EndDeepMIP(ii); 
        if any(isnan(depth))
            depth(1) = PETM.Start(ii);
            depth(2) = PETM.EndPeak(ii);
        end
    end
    
    depth = sort(depth);
    depthvar = PETM.DepthVariable(ii);
      
    if depthvar == "Stage"
        idx = find(Data.SiteName == PETM.SiteName(ii) & ...
               Data.ProxyType == PETM.ProxyType(ii) & ...
               strcmpi(Data.PublicationDOI, PETM.PublicationDOI(ii)) & ...
               Data.Stage == "Ypresian");
    else 
        idx = find(Data.SiteName == PETM.SiteName(ii) & ...
                   Data.ProxyType == PETM.ProxyType(ii) & ...
                   strcmpi(Data.PublicationDOI, PETM.PublicationDOI(ii)) & ...
                   Data.(depthvar) >= depth(1) & Data.(depthvar) <= depth(2));
    end
    
    if numel(idx) < 1
        fprintf('No data added: PETM site %s (%0.0f)\n',PETM.SiteName(ii),ii);
    else
        Data.Stage(idx) = "PETM";
    end

end

end
