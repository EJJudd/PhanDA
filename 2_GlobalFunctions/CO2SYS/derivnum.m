% derivnum()
% This subroutine computes partial derivatives of output carbonate variables 
% with respect to input variables (two), plus nutrients (two), ammonium and hydrogen sulfide,
% temperature and salinity, dissociation constants, and total boron.
%
% It uses central differences, introducing a small perturbation 
% (plus and minus of a delta) in one INPUT variable and computes the
% resulting induced change in each OUTPUT variable
%
% This subroutine has been modified from its original version (Orr et al.
% 2018) to allow for input variables of CO2, HCO3, and CO3, the inclusion
% of ammonium and hydrogen sulfide, and compatibility with CO2SYS.m(v3)
%
% After numerical tests, the small PERTURBATION (delta) is chosen
% as a fraction of a reference value as follows:
% * 1.e-3 (0.1%) for the equilibrium constants and total boron 
% * 1.e-6        for the pair of CO2 system input variables (PAR1, PAR2)
% * 1.e-4        for temperature and salinity
% * 1.e-3        for total dissolved inorganic P and Si (PO4, SI)
%
%**************************************************************************
%
%  **** SYNTAX:
%  [deriv, headers_der, units_der, headers_err, units_err]=...
%                  derivnum(VARID,PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%                           SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,...
%                           SI,PO4,NH4,H2S,pHSCALEIN,K1K2CONSTANTS,...
%                           KSO4CONSTANT,KFCONSTANT,BORON)
% 
%  **** SYNTAX EXAMPLES:
%  [der, headers, units] = derivnum('par1',2400,2200,1,2,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [der, headers, units] = derivnum('sit', 2400,   8,1,3,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [deriv, headers]      = derivnum(  'T',  500,   8,5,3,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [deriv]               = derivnum('S',2400,2000:10:2400,1,2,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [deriv]               = derivnum('K0',2400,2200,1,2,0:1:35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [deriv]               = derivnum('K1',2400,2200,1,2,35,0,25,0:100:4200,0,15,1,0,0,1,4,1,1,1)
%  
%**************************************************************************
%
% INPUT:
%
%   - VARID     :   character string to select the input variable for which
%                   derivatives will be taken with respect to. 
%                   This variable appears in denominator of each resulting derivative.
%                   = variable length, case insensitive, character code
%                     case 'par1'  :  Parameter 1 of the input pair (This is TAlk if PAR1TYPE is 1)
%                     case 'par2'  :  Parameter 2 of the input pair (This is TAlk if PAR2TYPE is 1)
%                     case 'sil', 'silt', 'tsil', 'silicate', or 'sit'     : Silicate concentration
%                     case 'phos', 'phost', 'tphos', 'phosphate', or 'pt'  : Phosphate concentration
%                     case 'amm', 'nh4', 'NH4'                    : Ammonia concentration
%                     case 'hyd', 'h2s', 'H2S'                    : Hydrogen Sulfide concentration
%                     case 't', 'temp' or 'temperature'           : temperature
%                     case 's', 'sal', or 'salinity'              : salinity
%                     case 'K0','K1','K2','Kb','Kw','Kspa', 'Kspc': dissociation constants 
%                     case 'bor': total boron
%
%   - all others :  same list of input parameters as in CO2SYS() subroutine (version 3, scalar or vectors)
%
%**************************************************************************%
%
% OUTPUT: 3 arrays	
%         a) an array containing the derivative of the following parameter (one row per sample):
%         b) the corresponding cell-array containing crudely formatted headers
%         c) the corresponding cell-array containing the units
%
%    POS  PARAMETER         UNIT
%
%    01 - TAlk              (umol/kgSW)
%    02 - TCO2              (umol/kgSW)
%    03 - [H+] in           (umol/kgSW)
%    04 - pCO2 in           (uatm)
%    05 - fCO2 in           (uatm)
%    06 - HCO3 in           (umol/kgSW)
%    07 - CO3 in            (umol/kgSW)
%    08 - CO2 in            (umol/kgSW)
%    09 - RF in             ()
%    10 - OmegaCa in        ()
%    11 - OmegaAr in        ()
%    12 - xCO2 in           (ppm)
%    13 - [H+] out          ()
%    14 - pCO2 out          (uatm)
%    15 - fCO2 out          (uatm)
%    16 - HCO3 out          (umol/kgSW)
%    17 - CO3 out           (umol/kgSW)
%    18 - CO2 out           (umol/kgSW)
%    19 - RF out            ()
%    20 - OmegaCa out       ()
%    21 - OmegaAr out       ()
%    22 - xCO2 out          (ppm)
%
% * 'in'  refers to INPUT  conditions (TEMPIN, PRESIN) as for CO2SYS
%   'out' refers to OUTPUT conditions (TEMPOUT, PRESOUT)
%
%
% CAUTION: derivnum.m is NOT necessarily designed to take partial derivatives
%          of input variables, only computed variables relative to input
%          variables. However, those partial derivatives of input variables
%          are kept in the OUT conditions to maintain consistency
%          of the order of output between the different input pairs.
%          Generally, we advise not to use these derivatives of input
%          variables. However, in some cases their results appear accurate.
%          Use them at your own risk.
%

function [derivatives, headers, units, headers_err, units_err] = ...
        derivnum (VARID,PAR1,PAR2,PAR1TYPE,PAR2TYPE, SAL,TEMPIN, ...
                  TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S, ...
                  pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT, ...
                  KFCONSTANT,BORON)

    % For computing derivative with respect to Ks, one has to call CO2sys with a perturbed K
    % Requested perturbation is passed through the following global variables
    global PertK    % Id of perturbed K
    global Perturb  % perturbation

    % Input conditioning
    % ------------------

    VARID = upper(VARID);
    % Determine lengths of input vectors
    veclengths=[length(PAR1) length(PAR2) length(PAR1TYPE)...
                length(PAR2TYPE) length(SAL) length(TEMPIN)...
                length(TEMPOUT) length(PRESIN) length(PRESOUT)...
                length(SI) length(PO4) length(NH4) length(H2S)...
                length(pHSCALEIN) length(K1K2CONSTANTS)...
                length(KSO4CONSTANT) length(KFCONSTANT)...
                length(BORON)];

    if length(unique(veclengths))>2
            disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
    end

    % Make column vectors of all input vectors
    PAR1         =PAR1         (:);
    PAR2         =PAR2         (:);
    PAR1TYPE     =PAR1TYPE     (:);
    PAR2TYPE     =PAR2TYPE     (:);
    SAL          =SAL          (:);
    TEMPIN       =TEMPIN       (:);
    TEMPOUT      =TEMPOUT      (:);
    PRESIN       =PRESIN       (:);
    PRESOUT      =PRESOUT      (:);
    SI           =SI           (:);
    PO4          =PO4          (:);
    NH4          =NH4          (:);
    H2S          =H2S          (:);
    pHSCALEIN    =pHSCALEIN    (:);
    K1K2CONSTANTS=K1K2CONSTANTS(:);
    KSO4CONSTANT =KSO4CONSTANT (:);
    KFCONSTANT   =KFCONSTANT   (:);
    BORON        =BORON        (:);

    % Find the longest column vector:
    ntps = max(veclengths);

    % Populate column vectors
    PAR1(1:ntps,1)          = PAR1(:)          ;
    PAR2(1:ntps,1)          = PAR2(:)          ;
    PAR1TYPE(1:ntps,1)      = PAR1TYPE(:)      ;
    PAR2TYPE(1:ntps,1)      = PAR2TYPE(:)      ;
    SAL(1:ntps,1)           = SAL(:)           ;
    TEMPIN(1:ntps,1)        = TEMPIN(:)        ;
    TEMPOUT(1:ntps,1)       = TEMPOUT(:)       ;
    PRESIN(1:ntps,1)        = PRESIN(:)        ;
    PRESOUT(1:ntps,1)       = PRESOUT(:)       ;
    SI(1:ntps,1)            = SI(:)            ;
    PO4(1:ntps,1)           = PO4(:)           ;
    SI(1:ntps,1)            = SI(:)            ;
    PO4(1:ntps,1)           = PO4(:)           ;
    NH4(1:ntps,1)           = NH4(:)           ;
    H2S(1:ntps,1)           = H2S(:)           ;
    pHSCALEIN(1:ntps,1)     = pHSCALEIN(:)     ;
    K1K2CONSTANTS(1:ntps,1) = K1K2CONSTANTS(:) ;
    KSO4CONSTANT(1:ntps,1)  = KSO4CONSTANT(:)  ;
    KFCONSTANT(1:ntps,1)    = KFCONSTANT(:)    ;
    BORON(1:ntps,1)         = BORON(:)         ;

    % BASELINE:
    % --------
    
    carb = CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);
    % Compute [H+] in µmol/KgSW
    if (ndims(carb) == 2)
        Hin = 10.^(-carb(:,3)) * 1.e6;
        Hout = 10.^(-carb(:,21)) * 1.e6;
        carb = horzcat(carb(:,1:2), Hin, carb(:,3:20), Hout, carb(:,21:end));
    else
        Hin = 10.^(-carb(3)) * 1.e6;
        Hout = 10.^(-carb(21)) * 1.e6;
        carb = horzcat(carb(1:2), Hin, carb(3:20), Hout, carb(21:end));
    end    

    % Compute two slightly different values for input
    % -----------------------------------------------

    % Default values : not perturbed
    % no change in PAR1 and PAR2
    PAR11 = PAR1; PAR12 = PAR1;
    PAR21 = PAR2; PAR22 = PAR2;
    % no change in Sil total and Phos total (except for d/dSi, d/dPO4)
    SI1 = SI;    SI2 = SI;
    PO41  = PO4; PO42  = PO4;
    % no change in Amm total and H2S total (except for d/dAmm, d/dH2S)
    NH41  = NH4; NH42 = NH4;
    H2S1  = H2S; H2S2  = H2S;
    % no change in T and S (except for d/dS and d/dT, respectively)
    TEMP1 = TEMPIN; TEMP2 = TEMPIN;
    SAL1 = SAL; SAL2 = SAL;
	% no change in TEMPOUT except for d/dT (see why further below)	
    TEMPOUT1 = TEMPOUT; TEMPOUT2 = TEMPOUT;

    % Create empty vector for abs_dx (absolute delta)
    abs_dx  = nan(ntps,1);

    % Flag for dissociation constant as perturbed variable 
    flag_dissoc_K = 0;   % False
    % names of 7 dissociation constants and variable for total boron 
    K_names = {'K0', 'K1', 'K2', 'KB', 'KW', 'KSPA', 'KSPC', 'BOR', 'CAL'};    

    % Units for derivatives and for errors
    units_at = {'umol';'umol';'nmol';'total scale';'uatm kg';'uatm kg';'umol';'umol';...
                 'umol';'kg';'kg';'kg';'ppm kg';...
                 'nmol';'uatm kg';'uatm kg';'umol';'umol';...
                 'umol';'kg';'kg';'kg';'ppm kg';
                };
 
    units_kg  = {'umol/kg';'umol/kg';'nmol/kg';'total scale';'uatm';'uatm';'umol/kg';'umol/kg';...
                 'umol/kg';' ';' ';' ';'ppm';...
                 'nmol/kg';'uatm';'uatm';'umol/kg';'umol/kg';...
                 'umol/kg';' ';' ';' ';'ppm';
                };

    units_k = units_kg;
 
    units_pco2 = units_kg;
    
    units_err = units_kg;
 
    switch VARID
        case K_names
            flag_dissoc_K = 1;
            % Approximate values for K0, K1, K2, Kb, Kspa, Kspc, and CAL
            % They will be used to compute an absolute perturbation
            % value on these constants
            K = [0.034, 1.2e-06, 8.3e-10, 2.1e-09, 6.1e-14, 6.7e-07, ...
                 4.3e-07, 0.0004157, 0.0103];
            % Choose value of absolute perturbation
            [is_in_K_names, index] = ismember(VARID, K_names);
            perturbation = K(index) * 1.e-3;   % 0.1 percent of Kx value
            abs_dx = 2 * perturbation;
            denom_headers = VARID;
            denom_units = ' ';
            units = units_k;
            
        case {'PAR1', 'VAR1'}      % PAR1 (first variable of input pair) is perturbed
            % Define a relative delta
            delta = 1.e-6;
          
            PAR1ref = nan(size(PAR1TYPE));
            units = units_kg;
            denom_headers = '<PAR1>';
            denom_units = '<PAR1 units>';
                t1=PAR1TYPE==1;
                  PAR1ref(t1) = 2300.; % umol/kg (global surface average, Orr et al., 2017)
                t2=PAR1TYPE==2;
                  PAR1ref(t2) = 2000.; % umol/kg (global surface average, Orr et al., 2017)
                t3=PAR1TYPE==3;
                  PAR1ref(t3) = 1.0e-8; % mol/kg (equivalent to pH=8.0)
                t4=PAR1TYPE==4;
                  PAR1ref(t4) = 400.;  % uatm
                t5=PAR1TYPE==5;
                  PAR1ref(t5) = 400.; % uatm
                t6=PAR1TYPE==6;
                  PAR1ref(t6) = 1790.; % umol/kg
                t7=PAR1TYPE==7;
                  PAR1ref(t7) = 200.; % umol/kg
                t8=PAR1TYPE==8;
                  PAR1ref(t8) = 10.; % umol/kg
      
           % cases where first input variable is pH
            F = (PAR1TYPE == 3);
            H = 10.^(-PAR1(F)) ; % [H+] in mol/kg
            % Change slightly [H+]
            H1 = H - PAR1ref(F)*delta;
            H2 = H + PAR1ref(F)*delta;
            PAR11(F) = -log10(H1) ;
            PAR12(F) = -log10(H2) ;
            abs_dx(F) = (H2 - H1) * 1e9; % now in nmol/kg
           
            G = ~F;
            % Change slightly PAR1
            PAR11(G) = PAR1(G) - PAR1ref(G)*delta;
            PAR12(G) = PAR1(G) + PAR1ref(G)*delta;
            abs_dx(G) = PAR12(G) - PAR11(G);

      case {'PAR2', 'VAR2'}    % PAR2 (second variable of input pair) is perturbed
            % Define a relative delta
            delta = 1.e-6;
            
            PAR2ref = nan(size(PAR2TYPE));
            units = units_kg;
            denom_headers = '<PAR2>';
            denom_units = '<PAR2 units>';
                t1=PAR2TYPE==1;
                  PAR2ref(t1) = 2300.; % umol/kg (global surface average, Orr et al., 2017)
                t2=PAR2TYPE==2;
                  PAR2ref(t2) = 2000.; % umol/kg (global surface average, Orr et al., 2017)
                t3=PAR2TYPE==3;
                  PAR2ref(t3) = 1.0e-8; % mol/kg (equivalent to pH=8.0)
                t4=PAR2TYPE==4;
                  PAR2ref(t4) = 400.;  % uatm
                t5=PAR2TYPE==5;
                  PAR2ref(t5) = 400.; % uatm
                t6=PAR2TYPE==6;
                  PAR2ref(t6) = 1790.; % umol/kg
                t7=PAR2TYPE==7;
                  PAR2ref(t7) = 200.; % umol/kg
                t8=PAR2TYPE==8;
                  PAR2ref(t8) = 10.; % umol/kg
            
            
            % cases where second input variable is pH
            F = (PAR2TYPE == 3);
            H = 10.^(-PAR2(F)) ; % H+ in mol/kg
            % Change slightly [H+]
            H1 = H - PAR2ref(F)*delta;
            H2 = H + PAR2ref(F)*delta;
            PAR21(F) = -log10(H1) ;
            PAR22(F) = -log10(H2) ;
            abs_dx(F) = (H2 - H1) * 1e9;

            G = ~F;
            % Change slightly PAR2
            PAR21(G) = PAR2(G) - PAR2ref(G)*delta;
            PAR22(G) = PAR2(G) + PAR2ref(G)*delta;
            abs_dx(G) = PAR22(G) - PAR21(G);

       case {'SIL', 'TSIL', 'SILT', 'SILICATE', 'SIT'}    % Sil total
            % Define a relative delta
            delta = 1.e-3;
            SIref = 7.5; % umol/kg (global surface average, Orr et al., 2018)
            % Change slightly SI
            SI1 = SI - SIref*delta;
            SI2 = SI + SIref*delta;
            abs_dx = SI2 - SI1;
            denom_headers = 'Sit';
            denom_units = 'umol';
            units = units_at;
       case {'PHOS', 'TPHOS', 'PHOST', 'PHOSPHATE', 'PT'}    % Phos total
            % Define a relative delta
            delta = 1.e-3;
            PO4ref = 0.5; % umol/kg (global surface average, Orr et al., 2018)
            % Change slightly PO4
            PO41 = PO4 - PO4ref*delta;
            PO42 = PO4 + PO4ref*delta;
            abs_dx = PO42 - PO41;
            denom_headers = 'Pt';
            denom_units = 'umol';
            units = units_at;
       case {'AMM','NH4'}    % Amm total
            % Define a relative delta
            delta = 1.e-3;
            NH4ref = 1;
            % Change slightly NH4
            NH41 = NH4 - NH4ref*delta;
            NH42 = NH4 + NH4ref*delta;
            abs_dx = NH42 - NH41;
            denom_headers = 'NH4t';
            denom_units = 'umol';
            units = units_at;
       case {'HYD','H2S'}    % H2S total
            % Define a relative delta
            delta = 1.e-3;
            H2Sref = 1;
            % Change slightly H2S
            H2S1 = H2S - H2Sref*delta;
            H2S2 = H2S + H2Sref*delta;
            abs_dx = H2S2 - H2S1;
            denom_headers = 'H2St';
            denom_units = 'umol';
            units = units_at;
       case {'T', 'TEMP', 'TEMPERATURE'}    % Temperature
            % Define a relative delta
            delta = 1.e-4;
            TEMPref = 18.; % global surface mean (C)
            % Change slightly temperature
            TEMP1 = TEMPIN - TEMPref*delta;
            TEMP2 = TEMPIN + TEMPref*delta;
            abs_dx = TEMP2 - TEMP1;
            % When computing d/dT	
            % TEMPOUT changes must be the same as TEMPIN changes, e.g, for	
            % derivnum 'OUT' & 'IN' results to be the same when TEMPOUT=TEMPIN	
            TEMPOUT1 = TEMPOUT - TEMPref*delta;	
            TEMPOUT2 = TEMPOUT + TEMPref*delta;
            denom_headers = 'T';
            denom_units = 'C';
            units = units_kg;
       case {'S', 'SAL', 'SALINITY'}    % Salinity
            % Define a relative delta
            delta = 1.e-4;
            SALref = 35.;
            % Change slightly salinity
            SAL1 = SAL - SALref*delta;
            SAL2 = SAL + SALref*delta;
            abs_dx = SAL2 - SAL1;
            denom_headers = 'S';
            denom_units = 'psu';
            units = units_kg;
    end
    
    % PERTURBATION:
    % -------------

    % Point 1: (one dissociation constant or PAR1, PAR2, T or S is somewhat smaller)
    % if perturbed variable is a dissociation constant
    if (flag_dissoc_K)
        PertK = upper(VARID);
        Perturb = -perturbation;
    end
    cdel1 = CO2SYS ( ...
        PAR11,PAR21,PAR1TYPE,PAR2TYPE,SAL1,TEMP1,TEMPOUT1,PRESIN,PRESOUT,SI1,PO41,NH41,H2S1,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);
    % Compute [H+]
    if (ndims(cdel1) == 2)
        Hin = 10.^(-cdel1(:,3))   * 1.e9; % to show H+ results in nmol/kg
        Hout = 10.^(-cdel1(:,21)) * 1.e9;
        cdel1 = horzcat(cdel1(:,1:2), Hin, cdel1(:,3:20), Hout, cdel1(:,21:end));
    else
        Hin = 10.^(-cdel1(3))     * 1.e9;
        Hout = 10.^(-cdel1(21))   * 1.e9;
        cdel1 = horzcat(cdel1(1:2), Hin, cdel1(3:20), Hout, cdel1(21:end));
    end

    % Point 2: (one dissociation constant or PAR1, PAR2, T or S is somewhat bigger)
    % if perturbed variable is a dissociation constant
    if (flag_dissoc_K)
        PertK = upper(VARID);
        Perturb = perturbation;
    end
    cdel2 = CO2SYS ( ...
        PAR12,PAR22,PAR1TYPE,PAR2TYPE,SAL2,TEMP2,TEMPOUT2,PRESIN,PRESOUT,SI2,PO42,NH42,H2S2,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);
    % Compute [H+]
    if (ndims(cdel2) == 2)
        % Computed variable H+ (does not affect other computed
        % variables, i.e., it is the numerator of the derivative)
        Hin = 10.^(-cdel2(:,3))   * 1.e9; % to show H+ results in nmol/kg 
        Hout = 10.^(-cdel2(:,21)) * 1.e9;
        cdel2 = horzcat(cdel2(:,1:2), Hin, cdel2(:,3:20), Hout, cdel2(:,21:end));
    else
        Hin = 10.^(-cdel2(3))     * 1.e9;
        Hout = 10.^(-cdel2(21))   * 1.e9;
        cdel2 = horzcat(cdel2(1:2), Hin, cdel2(3:20), Hout, cdel2(21:end));
    end

    % if perturbed variable is a dissociation constant
    if (flag_dissoc_K)
        PertK = ''; %% Return to normal
    end
    
    % Drop unnecessary columns
    % Note: columns (04 - pHin) and (20 - pHout) drop    
    % Keep only the following columns
    %    01 - TAlk              (umol/kgSW)
    %    02 - TCO2              (umol/kgSW)
    %    03 - [H+] in           (nmol/kgSW)  lastly added
    %    05 - pCO2 in           (uatm)
    %    06 - fCO2 in           (uatm)
    %    07 - HCO3 in           (umol/kgSW)
    %    08 - CO3 in            (umol/kgSW)
    %    09 - CO2 in            (umol/kgSW)
    %    17 - RF in             ()
    %    18 - OmegaCa in        ()
    %    19 - OmegaAr in        ()
    %    20 - xCO2 in           (ppm)
    %    22 - [H+] out          (nmol/kgSW)  lastly added
    %    24 - pCO2 out          (uatm)
    %    25 - fCO2 out          (uatm)
    %    26 - HCO3 out          (umol/kgSW)
    %    27 - CO3 out           (umol/kgSW)
    %    28 - CO2 out           (umol/kgSW)
    %    36 - RF out            ()
    %    37 - OmegaCa out       ()
    %    38 - OmegaAr out       ()
    %    39 - xCO2 out          (ppm)
    keep = [1 2 3 5 6 7 8 9 17 18 19 20 22 24 25 26 27 28 36 37 38 39];
    
    % We will drop also some column headers
    headers = {'TAlk';'TCO2';'Hin';'pHin';'pCO2in';'fCO2in';'HCO3in';'CO3in';...
        'CO2in';'RFin';'OmegaCAin';'OmegaARin';'xCO2in';...
        'Hout';'pCO2out';'fCO2out';'HCO3out';'CO3out';...
        'CO2out';'RFout';'OmegaCAout';'OmegaARout';'xCO2out'};
    headers_err = headers;
    %units = {'umol';'umol';'nmol';'total scale';'uatm';'uatm';'umol';'umol';...
    %    'umol';' ';' ';'ppm';...
    %    'nmol';'uatm';'uatm';'umol';'umol';...
    %    'umol';' ';' ';'ppm';
    %    };
    % Initially, keep all headers except 'pHin'
    keep_head =  [1:3 5:23];

%     **** This was previously uncommented to eliminate partial derivatives
%     **** of input variables
%
%     % if all parameter PAR1 are of the same type
%     if all(PAR1TYPE == PAR1TYPE(1))
%         % Determine column number of PAR1
%         if PAR1TYPE(1) <= 3 % TAlk, TCO2 or pH
%             % By design of CO2sys, PARTYPE is equal to column number
%             col_number = PAR1TYPE(1);
%         else
%             % Because there is an extra column: [H+]
%             col_number = PAR1TYPE(1) + 1;
%         end
%         % Exclude input parameters PAR1
%         A = (keep ~= col_number);
%         keep = keep (A);
%         A = (keep_head ~= col_number);
%         keep_head = keep_head (A);
%     end
%     % if all parameter PAR2 are of the same type
%     if all(PAR2TYPE == PAR2TYPE(1))
%         % Determine column number of PAR1
%         if PAR2TYPE(1) <= 3 % TAlk, TCO2 or pH
%             % By design of CO2sys, PARTYPE is equal to column number
%             col_number = PAR2TYPE(1);
%         else
%             % Because there is an extra column: [H+]
%             col_number = PAR2TYPE(1) + 1;
%         end
%         % Exclude input parameters PAR2
%         A = (keep ~= col_number);
%         keep = keep (A);
%         A = (keep_head ~= col_number);
%         keep_head = keep_head (A);
%     end
    
    cdel1 = cdel1(:,keep);
    cdel2 = cdel2(:,keep);
    
    headers = headers(keep_head);
    units = units(keep_head);
    
    headers_err = headers_err(keep_head);
    units_err = units_err(keep_head);
    
    % concatenate strings to give each header its proper form, a partial derivative
    headers = strcat('d', headers, '/', 'd', denom_headers);
    % concatenate in an analogous way for the units
    units   = strcat('(',units, '/', denom_units,')');
    
    % Centered difference
    dy = (cdel2 - cdel1);

    % Compute ratio dy/dx
    if (isscalar(abs_dx))
        derivatives = dy ./ abs_dx;
    else
        derivatives = bsxfun(@rdivide, dy, abs_dx);
    end
    
%     **** This was previously uncommented to mask derivatives of
%     **** parameters that are directly related to input variables
%
%     % Mask values that should not be used with NaN (e.g., dHout/dT when PAR1 or PAR2 is pH)
%     switch VARID
%         case {'T', 'TEMP', 'TEMPERATURE'}
%             % For PAR1TYPE or PAR2TYPE = 3 (pH is input) make dHout/dT value a NaN
%             F = (PAR1TYPE==3 | PAR2TYPE==3); % either CO2 system input variable is pH
%             [is_in_headers, idx] = ismember('dHout/dT', headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%             
%             % For PAR1TYPE or PAR2TYPE = 4 or 5 (pCO2 or fCO2 is input) make relevant d/dT values NaNs
%             F = (PAR1TYPE==4 | PAR2TYPE==4); %when pCO2 or fCO2 is input var
%             % masknan = {'dfCO2in/dT' 'dxCO2in/dT' 'dpCO2out/dT' 'dfCO2out/dT' 'dxCO2out/dT'};
%             masknan = {'dpCO2out/dT'};
%             [is_in_headers, idx] = ismember(masknan, headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%             
%             % For PAR1TYPE or PAR2TYPE = 5 (fCO2 is input) make relevant d/dT values NaNs
%             F = (PAR1TYPE==5 | PAR2TYPE==5); %when fCO2 is input var
%             % masknan = {'dfCO2in/dT' 'dxCO2in/dT' 'dpCO2out/dT' 'dfCO2out/dT' 'dxCO2out/dT'};
%             masknan = {'dfCO2out/dT'};
%             [is_in_headers, idx] = ismember(masknan, headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%         
%         case {'PAR1', 'VAR1'}
%             % For pH-pCO2 or pH-fCO2 pair (PAR1 is pH) 
%             F = (PAR1TYPE==3 & (PAR2TYPE==4 | PAR2TYPE==5));
%             masknan = {'dHout/dH' 'dpCO2out/dH' 'dfCO2out/dH'};
%             [is_in_headers, idx] = ismember(masknan, headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%         
%         case {'PAR2', 'VAR2'}
%             % For pCO2-pH or fCO2-pH pair (PAR2 is pH)
%             F = ((PAR1TYPE==4 | PAR1TYPE==5) & PAR2TYPE==3);
%             masknan = {'dHout/dH' 'dpCO2out/dH' 'dfCO2out/dH'};
%             [is_in_headers, idx] = ismember(masknan, headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%     end
    
end
