function [Yuse, Yeuse, Ruse, proxytype, paleolat, paleolon] = assembleYYeR( ...
    UPD, Y, Ye, Rvals, Assumptions, pHcorr, snowballcorr, Rmethod, a)


fn = fieldnames(Y);
Yuse = []; Yeuse = []; Ruse = []; 
proxytype = []; paleolat = []; paleolon = [];

for ii = 1:numel(fn)
    
    % Amass Y
    Yuse = [Yuse; Y.(fn{ii})];
    proxytype = [proxytype; repmat(string(fn(ii)), numel(Y.(fn{ii})), 1)];
    paleolat = [paleolat; UPD.(fn{ii}).PaleoLat];
    paleolon = [paleolon; UPD.(fn{ii}).PaleoLon];
    
    % Amass Ye
    Yetemp = Ye.(fn{ii});
    % Add in snowball correction if applicable
    if contains(fn(ii),{'d18a','d18c','d18p'})
        if snowballcorr
            Yetemp = Yetemp + Assumptions.global.d18Osw.Snowball;
        end
    end
    % Add in pH correction if applicable
    if contains(fn(ii),{'d18a','d18c'})
        if pHcorr == "ens"
            Yetemp = Yetemp + Assumptions.local.(fn{ii}).d18OpHens;
        elseif pHcorr == "rec"
            Yetemp = Yetemp + Assumptions.local.(fn{ii}).d18OpHrec;
        end
    end
    Yeuse = [Yeuse; Yetemp];

    % Amass R
    if Rmethod == "percentile"
        Rtemp = r_percentile(UPD, Rvals, fn{ii}, a);
    else
        Rtemp = r_fixed(UPD, Rvals, fn{ii}, a, Rmethod);
    end
    
    
    Ruse = [Ruse; Rtemp];
    
end

end