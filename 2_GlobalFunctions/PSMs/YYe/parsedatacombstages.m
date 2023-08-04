function [proxydata,uniqueproxydata] = parsedatacombstages( ...
    Data, Preferences, GTS, stageno)

stagename = GTS.Stage(stageno);

if any(stageno == Preferences.combstages(1,:)) 
    changestage = GTS.Stage(...
        Preferences.combstages(2,Preferences.combstages(1,:) == stageno));
elseif any(stageno == Preferences.combstages(2,:))
    changestage = GTS.Stage(...
        Preferences.combstages(1,Preferences.combstages(2,:) == stageno));
else
    changestage = '';
end

Data.Stage(string(Data.Stage) == changestage) = {stagename};
[proxydata,uniqueproxydata] = ...
    parsedata(Data, Preferences, stagename); 


end