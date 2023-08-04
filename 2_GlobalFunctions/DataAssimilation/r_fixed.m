function Rtemp = r_fixed(UPD, Rvals, proxytype, a, Rmethod)

Rtemp = NaN(numel(UPD.(proxytype).N),1);
load("GTS2020_PETM.mat")

% Select which d18acmacro & mg values to use
if  proxytype == "mg"
    if a <= find(GTS.Stage == "Langhian")
        rowuse = 1;
    else
        rowuse = 2;
    end
end

% Assign value
if Rmethod == "high"
    val = 3;
elseif Rmethod == "medium"
    val = 2;
elseif Rmethod == "low"
    val = 1;
end

% Assign Rtemp
if proxytype == "d18cforam"
    Rtemp(isnan(Rtemp)) = Rvals.d18cforam(val);
elseif proxytype == "d18cmacro" || proxytype == "d18a"
    Rtemp(isnan(Rtemp)) = Rvals.d18acmacro(val);
elseif proxytype == "d18p"
    Rtemp(isnan(Rtemp) & UPD.(proxytype).AnalyticalMethod == "IRMS") ...
        = Rvals.d18p(1,val);
    Rtemp(isnan(Rtemp) & UPD.(proxytype).AnalyticalMethod == "SIMS") ...
        = Rvals.d18p(2,val);
elseif proxytype == "mg"
    Rtemp(isnan(Rtemp)) = Rvals.mg(rowuse, val);
elseif proxytype == "tex"
    Rtemp(isnan(Rtemp)) = Rvals.tex(val);
elseif proxytype == "uk"
    Rtemp(isnan(Rtemp)) = Rvals.uk(val);
end
   


end