function d18Omin_out = temp2d18O(Tc, eq, mineralogy, d18Osw)

% Function to forward model d18O of a mineral (carbonate or phosphate), 
% given a specified temperature (celcius), equation, and mineralogy 
% (calcite, aragonite, or phosphate), and d18Osw

% Calcite and aragonite values are framed in the context of the
% fractionation factor (alpha), while the phosphate equations are based on
% raw transfer functions
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


% (1) convert temperature to kelvin
Tk = cel2kel(Tc);


% (2) Assign values A, B, & C to calculate fractionation factor, alpha
%   All calcite/aragonite alpha equations follow the same basic formula:
%     alpha = e ^ (A * B + C)
%   In most equations, 
%     B = Tk^-1
%   However in a a some equations (e.g., Grossman and Ku, 1986 and modified 
%   by Dettman et al., 1999), 
%     B = 1000*Tk^-2
%   Phosphate equations are based on the genernalized equation:
%     Tc = A + B * (d18min + C - d18sw)
%   where A & B are constants, and C is the correction to a NBS120c value
%   of 21.7 per mille. Rearranged in terms of d18min yields:
%     d18min = (Tc - A) / B - C + d18sw
%   For consistency, below, A is framed as Tc - A


if eq == "Daeron2019" && mineralogy == "calcite"
    A = 17.57;
    B = 1./Tk;
    C = -0.02913;
elseif eq == "Kim1997" && mineralogy == "calcite"
    A = 18.04;
    B = 1./Tk;
    C = -0.03218;
elseif eq == "Grossman1986" && mineralogy == "aragonite"
    A = 2.559;
    B = 1000./(Tk.^2);
    C = .000715;
elseif eq == "Kim2007" && mineralogy == "aragonite"
    A = 17.88;
    B = 1./Tk;
    C = -0.03114;
elseif eq == "Lecuyer2013" && mineralogy == "phosphate"
    A = Tc - 117.4;
    B = -4.50;
    C = 0;
elseif eq == "Puceat2010" && mineralogy == "phosphate"
    A = Tc - 118.7;
    B = -4.22;
    C = 22.6-21.7;
else
    warning('Invalid equation/minaralogy selected')
end

% (3) Solve for d18min
if mineralogy == "calcite" || mineralogy == "aragonite"
    % if calcite or aragonite, calculate alpha
    alpha_mw = exp(A .* B + C);
    % solve for d18Omin, following the equation:
    %   alpha_mw = (1000 + d18Omin) / (1000 + d18Osw)
    % where both d18O values are relative to VSMOW
    d18Omin_vsmow = alpha_mw .* (1000 + d18Osw) - 1000;
    % convert d18Omin from VSMOW to VPDB
    d18Omin_vpdb = vsmow2vpdb(d18Omin_vsmow);
    d18Omin_out = d18Omin_vpdb;
elseif mineralogy == "phosphate"
    % if phosphate, solve transfer function
    d18Omin_vsmow = A ./ B - C + d18Osw;
    d18Omin_out = d18Omin_vsmow;
end


end