function d18amacroEst = d18amacroPSM_forward(sst,d18sw)

% Grossman and Ku, 1986
d18amacroEst = (20.6-sst)./4.34-.27+d18sw;

end
