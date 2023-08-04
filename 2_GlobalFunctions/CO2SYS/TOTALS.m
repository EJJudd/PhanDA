function [Salts,Salts_Headers] = TOTALS(CO2SYSDATA,CO2SYSHEADERS)
%**************************************************************************
%
% This function calculates the total concentrations of conservative
% constituents from salinity using output from CO2SYS. Only total
% concentrations that are not output by the CO2SYS function are calculated.
%
% If HEADERS are provided, this function will work with output from any
% version of CO2SYS. If HEADERS are not provided, this function will work
% only with output from CO2SYS v3.0, found at:
% https://github.com/jonathansharp/CO2-System-Extd
%
% Citations:
% Relative mass fractions are taken from:
%    Millero, F.J., et al., 2008. Deep-Sea Research Part I 55, 50-72.
% Based on work from:
%    Riley, J.P. and Tongudai, M., 1967. Chemical Geology 2, 263-269.
%    Morris, A.W. and Riley, J.P., 1966. Deep-Sea Research 13, 669-705.
%    Culkin, F. and Cox, R.A., 1966. Deep-Sea Research 13, 789-804.
%    Carpenter, J.H. and Manella, M.E., 1973. J. Geophys. Res. 78, 3621-26.
%
% Function written by Jonathan Sharp
% University of South Florida
% June 16, 2020
%
%**************************************************************************

if exist('CO2SYSHEADERS','var')

   SALidx = find(strcmp(CO2SYSHEADERS,'SAL'));
    
% Cations (mol/kg-SW):
 
    % Sodium:
    Na = 0.5564924./22.989769.*(CO2SYSDATA(:,SALidx)./1.80655);
    % Potassium:
    K  = 0.0206000./39.0983.*(CO2SYSDATA(:,SALidx)./1.80655);
    % Magnesium:
    Mg = 0.0662600./24.305.*(CO2SYSDATA(:,SALidx)./1.80655);
    % Calcium:
    Ca = 0.0212700./40.087.*(CO2SYSDATA(:,SALidx)./1.80655);
    % Strontium:
    Sr = 0.0004100./87.62.*(CO2SYSDATA(:,SALidx)./1.80655);

% Anions (mol/kg-SW):

    % Chloride:
    Cl = 0.9989041./35.453.*(CO2SYSDATA(:,SALidx)./1.80655);
    % Bromide:
    Br = 0.0034730./79.904.*(CO2SYSDATA(:,SALidx)./1.80655);
   
else

% Cations (mol/kg-SW):
 
    % Sodium:
    Na = 0.5564924./22.989769.*(CO2SYSDATA(:,56)./1.80655);
    % Potassium:
    K  = 0.0206000./39.0983.*(CO2SYSDATA(:,56)./1.80655);
    % Magnesium:
    Mg = 0.0662600./24.305.*(CO2SYSDATA(:,56)./1.80655);
    % Calcium:
    Ca = 0.0212700./40.087.*(CO2SYSDATA(:,56)./1.80655);
    % Strontium:
    Sr = 0.0004100./87.62.*(CO2SYSDATA(:,56)./1.80655);

% Anions (mol/kg-SW):
 
    % Chloride:
    Cl = 0.9989041./35.453.*(CO2SYSDATA(:,56)./1.80655);
    % Bromide:
    Br = 0.0034730./79.904.*(CO2SYSDATA(:,56)./1.80655);
    
end

Salts = [Na,K,Mg,Ca,Sr,Cl,Br];
Salts_Headers = {'Na+','K+','Mg2+','Ca2+','Sr2+','Cl-','Br-'};

end
