function Rtemp = r_percentile(UPD, Rvals, proxytype, a)

load("GTS2020_PETM.mat")
% Determine the percentiles
Rpc = prctile(UPD.(proxytype).N,[100/3, 200/3]);
% Adjust in the event that the lower percentile is 1
Rpc(Rpc==1) = 1.01;

% Select which mg values to use
if proxytype == "mg"
    if a <= find(GTS.Stage == "Langhian")
        rowuse = 1;
    else
        rowuse = 2;
    end
end

% Initialize
Rtemp = NaN(numel(UPD.(proxytype).N),1);

if proxytype == "d18cforam"
    Rtemp(UPD.(proxytype).N < Rpc(1)) = Rvals.d18cforam(3);
    Rtemp(UPD.(proxytype).N > Rpc(2)) = Rvals.d18cforam(1);
    Rtemp(isnan(Rtemp)) = Rvals.d18cforam(2);
elseif proxytype == "d18cmacro" || proxytype == "d18a"
    Rtemp(UPD.(proxytype).N < Rpc(1)) = Rvals.d18acmacro(3);
    Rtemp(UPD.(proxytype).N > Rpc(2)) = Rvals.d18acmacro(1);
    Rtemp(isnan(Rtemp)) = Rvals.d18acmacro(2);    
elseif proxytype == "d18p"
    Rtemp(UPD.(proxytype).N < Rpc(1) & UPD.(proxytype).AnalyticalMethod == "IRMS") ...
        = Rvals.d18p(1,3);
    Rtemp(UPD.(proxytype).N < Rpc(1) & UPD.(proxytype).AnalyticalMethod == "SIMS") ...
        = Rvals.d18p(2,3);
    Rtemp(UPD.(proxytype).N > Rpc(2) & UPD.(proxytype).AnalyticalMethod == "IRMS") ...
        = Rvals.d18p(1,1);
    Rtemp(UPD.(proxytype).N > Rpc(2) & UPD.(proxytype).AnalyticalMethod == "SIMS") ...
        = Rvals.d18p(2,1);
    Rtemp(isnan(Rtemp) & UPD.(proxytype).AnalyticalMethod == "IRMS") ...
        = Rvals.d18p(1,2);
    Rtemp(isnan(Rtemp) & UPD.(proxytype).AnalyticalMethod == "SIMS") ...
        = Rvals.d18p(2,2);
elseif proxytype == "mg"
    Rtemp(UPD.(proxytype).N < Rpc(1)) = Rvals.mg(rowuse,3);
    Rtemp(UPD.(proxytype).N > Rpc(2)) = Rvals.mg(rowuse,1);
    Rtemp(isnan(Rtemp)) = Rvals.mg(2);
elseif proxytype == "tex"
    Rtemp(UPD.(proxytype).N < Rpc(1)) = Rvals.tex(3);
    Rtemp(UPD.(proxytype).N > Rpc(2)) = Rvals.tex(1);
    Rtemp(isnan(Rtemp)) = Rvals.tex(2);
elseif proxytype == "uk"
    Rtemp(UPD.(proxytype).N < Rpc(1)) = Rvals.uk(3);
    Rtemp(UPD.(proxytype).N > Rpc(2)) = Rvals.uk(1);
    Rtemp(isnan(Rtemp)) = Rvals.uk(2);
end
    
    
end