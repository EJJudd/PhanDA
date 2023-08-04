function AT = alphaT(T,S,pH)
%========================================================%
%
% File: alphaT.m
%
% Author:
%
% Richard E. Zeebe
% School of Ocean and Earth 
% Science and Technology 
% Department of Oceanography 
% University of Hawaii at Manoa
% 1000 Pope Road, MSB 504 
% Honolulu, HI 96822, USA
% http://www.soest.hawaii.edu/oceanography/faculty/zeebe.html
%
% Purpose:
%
% - MATLAB script to calculate oxygen isotope fractionation
%   between S = [CO2(aq.)] + [H2CO3] + [HCO3-] + [CO32-])
%   and H2O, described in:
%
%   Zeebe, R.E. "An expression for the overall oxygen 
%   isotope fractionation between the sum of dissolved 
%   inorganic carbon and water", Geochemistry, Geophysics, 
%   Geosystems, 2007.
%
% Modified by J. Tierney, 2021 into a function and vectorized.
%
% Disclaimer/Notes:
%
% - Individual fractionation factors between HCO3-, CO32-,
%   and H2O were measured between 15 and 40 deg C in NaHCO3
%   solutions by Beck et al., GCA, 69, 3493-3503, 2005. This 
%   should be kept in mind if alphaT is calculated at higher/
%   lower temperature or for seawater.
%
%
% - Carbonic acid dissociation constants themselves should be 
%   valid from S = 0-50 and TC = 0-50 deg C. They are based on 
%   the SEWATER-pH SCALE: Millero et al., Marine Chemistry, 100, 
%   80-94, 2006. pH input for S > 0 is hence on seawater scale. 
%   For scale conversion, see e.g. 
%
%   Zeebe and Wolf-Gladrow. "CO2 in Seawater: Equilibrium, 
%   Kinetics, Isotopes". pp. , Elsevier, 2001.
%
%   http://www.soest.hawaii.edu/oceanography/faculty/zeebe.html
%
% - Note that seawater alphaT (S=35) is not necessarily appropriate for 
%   understanding carbonate precipitation and 18O-fractionation
%   in marine organisms because the solution chemistry at the site of 
%   calcification is unlikely to resemble that of natural seawater.
%
%  
% - In this script, S is for salinity. The quantity
%   S = [CO2(aq.)] + [H2CO3] + [HCO3-] + [CO32-])
%   (calligraphic symbol S in the paper) is not used. 
%
%
% Check values for alphaT:
% 
% 1.03278391897034 (T = 25 deg C, S = 00, pH = 7)
% 1.03186728109924 (T = 25 deg C, S = 35, pH = 7)
% 1.04337976403670 (T = 15 deg C, S = 35, pH = 0)
% 1.02190868921546 (T = 40 deg C, S = 35, pH =14)
%
% 
% Updates:
% 
% 04/15/07 new
% 06/05/07 comments added
% 05/20/21 modified by J. Tierney into a function and vectorized.
%============================================%
%
% Input Section. Specify:
%
% 1. Temperature
% 2. Salinity
% 3. pH 
%
% and run. 
%
% ranges:
% T  = 15-40 deg C
% S   = 0-?
% pH  = 0-14


%============================================%
%
% Code below this line is for calculating alphaT
%
%============================================%
%ensure column vectors
T = T(:);
S = S(:);
pH = pH(:);

%---------------------------------------------
%
% Calculate CO2 constants. Millero et al. 2006
%
%---------------------------------------------

Tk = 273.15;
T  = T+Tk;

% K1
A1   =   13.41910.*sqrt(S) + 0.0331.*S - 5.33e-5.*S.*S;
B1   = -530.12300.*sqrt(S) - 6.1030.*S;
C1   =   -2.06950.*sqrt(S);

pK10 = -126.34048 + 6320.813./T + 19.568224*log(T);
pK1  = pK10 + A1 + B1./T + C1.*log(T);
K1   = 10.^(-pK1);

% K2
A2   =   21.0894.*sqrt(S) +  0.1248.*S - 3.687e-4.*S.*S;
B2   = -772.4830.*sqrt(S) - 20.0510.*S;
C2   = -  3.3336.*sqrt(S);

pK20 = -90.18333 + 5143.692./T + 14.613358*log(T);
pK2  = pK20 + A2 + B2./T + C2.*log(T);
K2   = 10.^(-pK2);

K1p  = 10.^(-2.98);


%---------------------------------------------
%
% Calculate concentrations of carbonate species
%
%---------------------------------------------

% Notation:
%
% co2  = [CO2(aq.)]
% co2s = [CO2(aq.)] + [H2CO3]

dic = 2.0e-3; % mol/kg (value irrelevant for alpha T) 
              % only used here for consistency with equations
              % in Zeebe & Wolf-Gladrow, 2001.

h   = 10.^(-pH);


co2s   = dic./(1+K1./h+K1.*K2./h./h);
co2    = co2s/(K1p+1);
h2co3  = K1p*co2;
co2    = co2s-h2co3;
hco3   = dic./(1+h./K1+K2./h);
co3    = dic./(1+h./K2+h.*h./K1./K2);

% short:
d = co2;   % d: dioxide (aq.)
a = h2co3; % a: acid
b = hco3;  % b: bicarbonate
c = co3;   % c: carbonate

V = 2.*d+3.*(a+b+c);

% z-ratios:
zd =   2.*d./V;
za =   3.*a./V;
zb =   3.*b./V;
zc =   3.*c./V;


%---------------------------------------------
%
% Isotopes. Beck et al., 2005
%
%---------------------------------------------

lnad1000 = 2.52e6./T./T+12.12; % alpha(CO2aq-H2O)
lnab1000 = 2.59e6./T./T+ 1.89; % alpha( HCO3-H2O)
lnac1000 = 2.39e6./T./T- 2.70; % alpha(  CO3-H20)
aa       = 1.0400; % alpha(H2CO3-H2O) ~no effect on alphaT
ad       = exp(lnad1000/1000);
ab       = exp(lnab1000/1000);
ac       = exp(lnac1000/1000);


RH2O = 2.e-3; % 2.e-3 value irrelevant over
              % several orders of magnitudes
 
Rd = ad.*RH2O;
Ra = aa.*RH2O;
Rb = ab.*RH2O;
Rc = ac.*RH2O;

rd = Rd./(1+Rd);
ra = Ra./(1+Ra);
rb = Rb./(1+Rb);
rc = Rc./(1+Rc);

% R_Total
RT =   (rd.*zd+ra.*za+rb.*zb+rc.*zc)./ ...
    (1-(rd.*zd+ra.*za+rb.*zb+rc.*zc));
  
% alpha_Total
AT   = RT./RH2O;

end