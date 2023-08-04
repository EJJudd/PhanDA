function d18pEst = d18pPSM_forward(sst,d18sw)

%Lécuyer et al., 2013
d18pEst = (117.4-sst)./4.50+d18sw;


% Kolodny et al., 1983
% d18pEst = (113.3-sst)./4.38+d18sw;

end
