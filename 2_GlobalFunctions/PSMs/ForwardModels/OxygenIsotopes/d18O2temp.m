function Tc = d18O2temp(d18Omin_vpdb, eq, mineralogy, d18Osw)

% Function to estimate temperature of precipitation of a mineral 
% (carbonate or phosphate), given a specified isotopic value (VPDB), 
% equation, and mineralogy (calcite, aragonite, or phosphate), and d18Osw
%
% Options for calcite mineralogy:
%   (a) Daeron2019 - Daeron et al., 2019 
%       (https://doi.org/10.1038/s41467-019-08336-5)
%       Best applied to slow growing calcite & belemnites
%   (b) Kim1997 - Kim & O'Niel, 1997
%       (https://doi.org/10.1016/S0016-7037(97)00169-5)
%       Recalculated for IUPAC AFF, following Bajnai & Tödter, 2021
%       (https://davidbajnai.github.io/isogeochem/)
%       Derived using synthetic calcite grown at temps between 0-40 C
% Options for aragonite mineralogy:
%   (a) Grossman1986 - Grossman & Ku, 1986
%       (https://doi.org/10.1016/0168-9622(86)90057-6)
%       Modified for VSMOW by Dettman et al., 1999
%       (https://doi.org/10.1016/S0016-7037(99)00020-4)
%   (b) Kim2007 - Kim et al., 2007
%       (https://doi.org/10.1016/j.chemgeo.2007.08.005)
% Options for phosphate mineralogy:
%   (a) Lecuyer2013 - Lécuyer et al., 2013
%       (https://doi.org/10.1016/j.chemgeo.2013.03.008)
%   (b) Puceat2010 - Pucéat et al., 2010
%       (https://doi.org/10.1016/j.epsl.2010.07.034)


% (1) convert carbonate d18Omin from PDB to VSMOW
if mineralogy ~= "phosphate"
    d18Omin_vsmow = vpdb2vsmow(d18Omin_vpdb);
else
    d18Omin_vsmow = d18Omin_vpdb;
end

% (2) calculate carbonate alpha term, following the equation:
%       alpha_mw = (1000+d18Omin)/(1000+d18Osw)
%     where both d18O values are relative to VSMOW
alpha_mw = (1000 + d18Omin_vsmow) ./ (1000 + d18Osw);

% (3) Assign values A, B, & C to calculate fractionation factor, alpha
%   All alpha equations follow the same basic formula:
%     alpha = e ^ (A * B + C)
%   In most equations, 
%     B = Tk^-1
%   However in a a some equations (e.g., Grossman and Ku, 1986 and modified 
%   by Dettman et al., 1999), 
%     B = 1000*Tk^-2
%   Phosphate equations are based on the genernalized equation:
%     Tc = A + B * (d18min + C - d18sw)
%   where A & B are constants, and C is the correction to a NBS120c value
%   of 21.7 per mille.

if strcmpi(eq, 'Daeron2019') && strcmpi(mineralogy, 'calcite')
    A = 17.57;
    C = -0.02913;
    b = @(x) 1./x;
elseif strcmpi(eq, 'Kim1997') && strcmpi(mineralogy, 'calcite')
    A = 18.04;
    C = -0.03218;
    b = @(x) 1./x;
elseif strcmpi(eq, 'Grossman1986') && strcmpi(mineralogy, 'aragonite')
    A = 2.559;
    C = .000715;
    b = @(x) sqrt(1000./x);
elseif strcmpi(eq, 'Kim2007') && strcmpi(mineralogy, 'aragonite')
    A = 17.88;
    C = -0.03114;
    b = @(x) 1./x;
elseif strcmpi(eq, 'Lecuyer2013') && strcmpi(mineralogy, 'phosphate')
    A = 117.4;
    C = 0;
    B = -4.5;
elseif strcmpi(eq, 'Puceat2010') && strcmpi(mineralogy, 'phosphate')
    A = 118.7;
    C = 22.6-21.7;
    B = -4.22;
else
    warning('Invalid equation/minaralogy selected')
end

% (4) Calculate temperature 
if mineralogy ~= "phosphate"
%   (a) Solve for B for carbonate data
    B = ( log(alpha_mw) - C ) / A;
%   (b) Solve for temperature, in Kelvin
    Tk = feval(b,B);   
%   (c) Convert temp from Kelvin to Celcius
    Tc = kel2cel(Tk);
else
%   (a) Solve linear equation in C
    Tc = A + B * (d18Omin_vsmow + C - d18Osw);
end



end